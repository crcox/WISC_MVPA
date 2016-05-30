function [voxDim,minXYZ] = voxelspacing(coordinates,varargin)
% voxelspacing  Find largest voxel dimensions that do not cause overlap
%   [voxDim,minXYZ] = voxelspacing(coordinates);
%   [voxDim,minXYZ] = voxelspacing(coordinates, varargin);
%
%   INPUT
%   coordinates : A cell array of coordinate spaces.
%   OPTIONS:
%   voxdim      : Initial guess of the voxel dimensions.
%   stepsize    : Amount by which to adjust voxel dimensions on each
%                 iteration.
%
%   OUTPUT
%   voxDim : The discovered voxel dimensions.
%   minXYZ : The smallest XYZ values across all coordinate spaces.

    p = inputParser();
    addRequired(p, 'coordinates');
    addOptional(p, 'voxdim', [], @isrow);
    addOptional(p, 'stepsize', 0.25, @isscalar);
    parse(p, coordinates, varargin{:});

    if iscell(p.Results.coordinates)
      XYZ = p.Results.coordinates;
    else
      XYZ = {p.Results.coordinates};
    end
    voxDim = p.Results.voxdim;
    STEP_SIZE = p.Results.stepsize;

    if isempty(voxDim)
        voxDim = ones(1,size(XYZ{1},2));
    end

    minXYZ = min(cell2mat(XYZ));
    XYZm = cellfun(@(x) bsxfun(@minus,x,min(x)), XYZ, 'unif', 0);
    N = cellfun(@(x) size(x,1), XYZ);
    n = numel(XYZ);

    ndim = length(voxDim);
    d = -1;
    breakloop = false;
    warning('off', 'DefineCommonGrid:roundto:duplicate');
    while true
        for i = 1:n
            xyz = XYZm{i};
            xyz = roundto(xyz, voxDim);
            xyz_u = unique(xyz,'rows');

            if size(xyz_u, 1) < N(i)
                breakloop = true;
                voxDim(dd) = voxDim(dd) - STEP_SIZE;
                break
            else
                dd = mod(d,ndim)+1;
                d = d + 1;
                voxDim(dd) = voxDim(dd) + STEP_SIZE;
            end
        end
        if breakloop
            break
        end
    end
    warning('on', 'DefineCommonGrid:roundto:duplicate');
end
