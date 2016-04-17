function [rowfilter, colfilter] = composeFilters(FILTERS,labels)

  labels = ascell(labels);

  % metadata.filter points to a structured array of filters.
  % First, force filters to a common orientation.
  for ii = 1:numel(FILTERS)
    FILTERS(ii).filter = forceRowVec(FILTERS(ii).filter);
  end
  
  % If necessary, combine filters that are mean to be 'or-ed'
  z = cellfun(@iscell, labels);
  if any(z)
      prototype = [...
          forceRowVec(fieldnames(FILTERS(1)));...
          cell(1,numel(fieldnames(FILTERS(1))))...
      ];
      N = numel(FILTERS);
      FILTERS(end+nnz(z)) = struct(prototype{:});
      to_or = labels(z);
      label_comb_list = cell(1,nnz(z));
      for i = 1:nnz(z);
          label_set = to_or{i};
          n = numel(label_set);
          label_set_tagged = cell(1,n+1);
          label_set_tagged(1:n) = label_set;
          label_set_tagged{end} = 'autogen';
          label_comb = strjoin(label_set_tagged, '_');
          label_comb_list{i} = label_comb;
          filter_selection = FILTERS(ismember({FILTERS(1:N).label}, label_set));
          filter_comb = any(cat(1, filter_selection.filter),1);
          FILTERS(N+i).label = label_comb;
          FILTERS(N+i).dimension = filter_selection(1).dimension;
          FILTERS(N+i).filter = filter_comb;
      end
      labels(z) = label_comb_list;
  end

  % Then select the filters
  z = ismember({FILTERS.label}, labels);

  filters_row = FILTERS(z & [FILTERS.dimension]==1);
  filters_col = FILTERS(z & [FILTERS.dimension]==2);
  rowfilter = all(cat(1, filters_row.filter),1);
  colfilter = all(cat(1, filters_col.filter),1);
end

function C = ascell(X)
  if ~iscell(X)
    C = {X};
  else
    C = X;
  end
end