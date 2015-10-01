The CMU IDM format contains a structure named "meta" which keeps the information
about a dataset that is not directly related to the experimental design. The format
assumes the data is a matrix with <n> examples (or time points) by <m> voxels.
and generally only a subset of the voxels in the volume is present. Hence, it's
necessary to keep a mapping between columns of the data matrix and position in 3D.
This is done with the following fields of meta:

dimx,dimy,dimz - the dimensions of the imaging volume

colToCoord - m x 3 matrix 

	   vector = coordToCol(i,:) places the 3D coordinates of voxel i in vector

coordToCol - <dimx> x <dimy> x <dimz> matrix

	   i = colToCoord(x,y,z) gives the column of voxel with 3D coordinates x, y and z
	   If i = 0 then the voxel is not in the data matrix.

indicesIn3D - m element vector

	   these are just the linear indices into a <dimx> x <dimy> x <dimz> matrix
	   of the voxels in the data matrix, i.e.

	   volume(meta.indicesIn3D) should return 1 through m

	   This is useful in order to load a vector of m elements
	   (e.g. the p-values of some test over each data matrix column
	   into a volume, for instance


	   volume = repmat(NaN,[dimx,dimy,dimz]); % creates a 3D volume with NaN values
	   volume(meta.indicesIn3D) = vector;     % places vector in that 3D volume


The easiest way to create the structure is, assuming you have a 3D matrix with a binary mask
of the voxels you are interested in, "mask" (can be created easily by reading in
an AFNI mask with afni_matlab, say):


meta.dimx = dimx; meta.dimy = dimy; meta.dimz = dimz;
meta.dimensions = [dimx dimy dimz];

meta.indicesIn3D = find( mask(:) );
m = length(meta.indiceIn3D);
[cx,cy,cz] = ind2sub(meta.dimensions,indicesIn3D);
meta.colToCoord = [cx,cy,cz];
meta.coordToCol = zeros(meta.dimensions);
meta.coordToCol(meta.indicesIn3D) = 1:m;


This can then be used to produce the two remaining fields

numberOfNeighbours - m x 1

	numberOfNeighbours(i) is the number of neighbours a voxel has

voxelsToNeighbours - m x (2*radius+1)^3

	neighbours = voxelsToNeighbours(i,1:numberOfNeighbours(i));

	are the neighbours of voxel i (not including i).
	Note that the row <i> in voxelsToNeighbours will be filled with numbers only
	if voxel <i> has neighbours in every possible direction. If not, the entries
	above numberOfNeighbours(i) are just NaN.

These two fields are produced using a function that takes the coordinate information
and a specific radius to consider (e.g. radius 1 means a voxel and its immediately
adjacent neighbours).

[meta.voxelsToNeighbours,meta.numberOfNeighbours] = neighboursWithinRadius(meta.colToCoord,1);

Note that there are various kinds of neighbourhood definitions. The code defaults to 'cubic',
where radius 1 includes all 26 voxels adjacent, but one could also specify 'spheric' to
get the definition in [Kriegeskorte et al 2006], e.g.

[meta.voxelsToNeighbours,meta.numberOfNeighbours] = neighboursWithinRadius(meta.colToCoord,1,'neighbourhoodType','spheric');

(WARNING: 'spheric' not implemented yet)
