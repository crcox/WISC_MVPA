function [windows,window,corners] = ndmovingwindow(dims, varargin)
%sub2ind([10,10,10],[1,2,1,2,1,2,1,2],[1,1,2,2,1,1,2,2],[1,1,1,1,2,2,2,2])
  p = inputParser;
  addRequired(p, 'dims');
  addParameter(p, 'size', 1, @isnumeric);
  addParameter(p, 'overlap', 0, @isnumeric);
  parse(p, dims, varargin{:});
  
  WindowSize = p.Results.size;
  Overlap = p.Results.overlap;
  dims = p.Results.dims;
  
  if nargin > 1
    assert(numel(dims) == numel(WindowSize))
  end
  
  extent = arrayfun(@(x) 1:x, WindowSize, 'Unif', 0);
  E = cell(1, numel(extent));
  [E{:}] = ndgrid(extent{:});
  E = cellfun(@(x) x(:), E, 'Unif', 0);
  
  if numel(dims) > 1
    window = sub2ind(dims, E{:})-1;
  else
    window = (1:WindowSize)';
  end
  
  corners = cell(numel(dims),1);
  d = [0, cumprod(dims(1:end-1))];
  for i = 1:numel(dims)
    corners{i} = (1:(WindowSize(i)-Overlap(i)):dims(i)) + d(i);
  end
  corners = cell2mat(corners);
  windows = bsxfun(@plus, window, corners(:)');
end
  
  