The CMU IDM format contains a structure named "meta" which keeps the information
about a dataset that is not directly related to the experimental design. The format
assumes the data is a matrix with <n> examples (or time points) by <m> voxels.
and generally only a subset of the voxels in 3D the volume is present. Hence, it's
necessary to keep a mapping between columns of the data matrix and position in 3D.
All the code in the toolbox assumes you are are working with data in this format.

In order to describe how to convert to this format, I'll assume you have
the following two MATLAB matrices

a) mask - a 3D (dimx x dimy x dimz) binary mask with 1 where brain voxels are and 0 everywhere else
b) data - a 4D (dimx x dimy x dimz x dimt) dataset matrix, with <dimt> 3D volumes

The steps are:


1) create a "meta" structure from the mask

[meta] = createMetaFromMask( mask );

more details on the contents of "meta" are below.
The part that will be used in the conversion is the field

meta.indicesIn3D

which contains the indices of the 1 voxels in the mask when the 3D volume is vectorized.


2) turn each volume in data into one example

[dimx,dimy,dimz,dimt] = size(data);
nvoxels = length(meta.indicesIn3D);

examples = zeros(dimt,nvoxels);

for t = 1:dimt

    volume = data(:,:,:,t);
    examples(t,:) = volume(meta.indicesIn3D);

end 

sometimes I have to transform several volumes into a single example (e.g. all the images in a block).


3) create "labels" and "labelsGroup"

This depends entirely on how you want to label examples relative to TRs, the main thing
is that it has to have as many entries as there are examples.

Many of the functions ask for "group labels". This is a set of labels that groups examples for
cross-validation, e.g. if you want to leave examples from one run out, this would simply be
the number of the run for each example. If you are doing things at a finer temporal grain,
a group could be a number of consecutive images that you don't want ending up split across
the train and test sets. Much as the "labels", this needs to have as many entries as there are
examples.

The functions that use this will likely perform some sort of cross-validation inside, e.g.
when creating a searchlight accuracy map.

and you are all set!

You should probably save these rather than recompute them on the fly, as building meta
takes a little while...


----------------------------------------------------------------------------------------------

If you are more curious about "meta", it is a structure with the following fields:

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


The way createMetaFromMask works, assuming you have a 3D matrix with a binary mask
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
