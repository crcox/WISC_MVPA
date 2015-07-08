function c = ndcoord(varargin)
    g = cell(1,length(varargin));
    [g{:}] = ndgrid(varargin{:});
    c = cell2mat(cellfun(@(x) x(:), g, 'Unif', 0));
end