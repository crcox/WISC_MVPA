function [rowfilter, colfilter] = composeFilters(FILTERS,labels,varargin)
  p = inputParser;
  p.KeepUnmatched = false;
  addParameter(p , 'logic', @all, @isfunction_handle);
  parse(p, varargin{:});
  combinelogic = p.Results.logic;

  labels = ascell(labels);

  % metadata.filter points to a structured array of filters.
  % First, force filters to a common orientation.
  for ii = 1:numel(FILTERS)
    FILTERS(ii).filter = forceRowVec(FILTERS(ii).filter);
  end

  % Then select the filters
  z = false(1,numel(FILTERS));
  for f = labels;
    z(strcmp(f, {FILTERS.label})) = true;
  end

  filters_row = FILTERS(z & [FILTERS.dimension]==1);
  filters_col = FILTERS(z & [FILTERS.dimension]==2);
  rowfilter = combinelogic(cat(1, filters_row.filter),1);
  colfilter = combinelogic(cat(1, filters_col.filter),1);
end
