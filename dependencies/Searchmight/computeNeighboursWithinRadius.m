%
% Find the neighbours of a voxel within a specified radius (cube) 
%
% In:
% - a #voxels x 3 matrix, where each row has the 3D coordinates of each voxel
% - radius
%
% Out:
% - a matrix with the indices of the neighbours, e.g.
%     row <v> contains a list with the <k> neighbours of voxel <v>
%     followed by 0s, up to <radius+2>^3 (27 for radius 1).
%
% - a vector with the count of neighbours, e.g.
%     entry <v> contains <k>, the # of neighbours in row <v> in the previous matrix
%
% Examples:
%
% History:
% - 2009 Feb 5 - fpereira@princeton.edu - created from previous code
%

function [voxelsToNeighbours,numberOfNeighbours] = computeNeighboursWithinRadius( varargin ) 

%
% Process parameters and figure out a few things
%

this = 'computeNeighboursWithinRadius';
if nargin >=2
  colToCoord = varargin{1};

  nVoxels = size(colToCoord,1);
  radius  = varargin{2};
  neighbourhoodType = 'cubic';
  
  idx = 3;
  while idx <= nargin
    argval = varargin{idx};
    switch argval
     case {'neighbourhoodType'}
      neighbourhoodType = varargin{idx}; idx = idx + 1;
    end
    idx = idx + 1;
  end

  % check things
  radius = round(radius);
  if radius < 1; fprintf('error: radius must be at least 1 or a larger integer\n');return; end
  switch neighbourhoodType
   case {'cubic'}
    % OK
   otherwise
    fprintf('error: neighbourhood type %s is not supported\n',neighbourhoodType);return;
  end
else
  eval(sprintf('help %s',this)); return;
end

% find the neighbours of each voxel
% - numberOfNeighbours: #voxels x 1 (number of neighbours for voxel in columnsToNeighbours)
% - voxelsToNeighbours: #voxels x radius^3
%   Only the first #ofneighbours positions are filled, e.g. the
%   neighbours of voxel <v> are:
%   voxelsToNeighbours(v,1:numberOfNeighbours(v));

[voxelsToNeighbours,numberOfNeighbours] = findColumnNeighboursRadius(colToCoord,radius,neighbourhoodType);
  
  
  
% Find the columns that contain neighbours of each voxel within a radius of r voxels
%
% Outputs:
% - columnsToNeighbours: #voxels x 26 (only the first #ofneighbours positions are filled)
% - # of neighbours: #voxels x 1 (number of neighbours for voxel in columnsToNeighbours)
%

function [columnsToNeighbours,numberOfNeighbours] = findColumnNeighboursRadius(colToCoord,radius,neighbourhoodType)
  nvoxels     = size(colToCoord,1);  
  nvoxelsSide = 2*radius + 1;
  nvoxelsCube = nvoxelsSide^3;
  
  DEBUG = 0;
  
  columnsToNeighbours = zeros(nvoxels,nvoxelsCube-1);
  numberOfNeighbours  = zeros(nvoxels,1);
  
  if radius == 0; return; end

  fprintf('findColumnNeighboursRadius: finding neighbours within radius %d\n',radius);
  
  for v = 1:nvoxels

    if ~rem(v,100); fprintf('%d ',v); end
    
    % compute the coordinatewise distances between the voxel and all others
    vcoords   = colToCoord(v,:);
    tmp       = repmat(vcoords,nvoxels,1);
    distances = abs(colToCoord - tmp);
    
    % now find the neighbours (voxels within the radius)
    % these are voxels that are within <radius> for all three coordinates
    mask = (distances <= radius);
    mask(v,:)   = [0 0 0]; % exclude the voxel itself
    tmp         = sum(mask,2);
    isneighbour = (tmp == 3);
    neighbours  = find(isneighbour);
    
    numberOfNeighbours(v) = length(neighbours);
    columnsToNeighbours(v,1:numberOfNeighbours(v)) = neighbours';
  
    if DEBUG
      if numberOfNeighbours(v)
	fprintf('%d with coords %d,%d,%d has neighbours:\n',v,vcoords);
	
	for n = 1:numberOfNeighbours(v)
	  neighbour       = columnsToNeighbours(v,n);
	  neighbourCoords = colToCoord(neighbour,:);
	  
	  fprintf('\t%d\t%d,%d,%d\n',neighbour,neighbourCoords);	
	end
      else
	fprintf('%d with coords %d,%d,%d has no neighbours\n',v,vcoords);
      end
    end
  end
  
  fprintf('\nfindColumnNeighboursRadius: finding neighbours done\n');
