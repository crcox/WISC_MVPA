function [rowfilter, colfilter] = composeFilters(filterset,labels,varargin)
  p = inputParser;
  p.KeepUnmatched = false;
  addParameter(p , 'logic', @all, @isfunction_handle);
  parse(p, varargin{:});
  combinelogic = p.Results.logic;

  labels = ascell(labels);

  % metadata.filter points to a structured array of filters.
  % First, force filters to a common orientation.
  for ii = 1:numel(filterset)
    filterset(ii).filter = forceRowVec(filterset(ii).filter);
  end

  % Then select the filters
  z = false(1,numel(filterset));
  for f = labels;
    z(strcmp(f, {filterset.label})) = true;
  end

  filters.row = filterset(z & [filterset.dimension]==1);
  filters.col = filterset(z & [filterset.dimension]==2);
  rowfilter = combinelogic(cat(1, filters.row.filter),1);
  colfilter = combinelogic(cat(1, filters.col.filter),1);
end
