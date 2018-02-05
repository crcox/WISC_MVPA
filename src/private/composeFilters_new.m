function [rowfilter, colfilter] = composeFilters_new(FILTERS,labels)
% COMPOSEFILTERS Generate a single row- and column-filter from multiple.
%
% FILTERS is an array of substructures of the standard METADATA structure. It
% has three fields:
%  1. label     - a unique string identifier for referencing the filter.
%  2. dimension - Indicates whether the filter is applied to rows or columns of
%                 the data matrix. Acceptable values are 1 (applies to rows) or
%                 2 (applies to column).
%  3. filter    - The filter itself, represented as a logical vector. Row
%                 filters must have as many elements as there are rows in the
%                 data matrix to which it will be applied, and vice versa for
%                 column filters. By extension, within a particular FILTERS
%                 structure, all filters of a given dimension must be of the
%                 same length.
%                 It does not matter whether filters are stored as row or
%                 column vectors. It is the responsibility of functions that
%                 reference this structure to check the orientation and
%                 transpose as necessary. This function will force all filters
%                 to be rowvectors, so that all() and any() can be applied by
%                 column (their default mode of operation).
%
% labels is a cell array of strings, corresponding to labels in the FILTERS
% structured array. As an example:
%
%   FILTERS(1) = struct('label','A','dimension',1,'filter',[1,1,1,0]);
%   FILTERS(2) = struct('label','B','dimension',1,'filter',[1,1,0,1]);
%   FILTERS(3) = struct('label','C','dimension',2,'filter',[1,1,1,0]);
%   FILTERS(4) = struct('label','D','dimension',2,'filter',[1,1,0,1]);
%   FILTERS(5) = struct('label','E','dimension',2,'filter',[1,0,1,1]);
%   FILTERS(6) = struct('label','F','dimension',1,'filter',[0,0,1,1]);
%
%   labels = {'A','B',{'C','D'},'E'};
%   [rowfilter,colfilter] = composeFilters(FILTERS, labels)
%
%   rowfilter =
%     [1,1,0,0]
%
%   colfilter =
%     [1,0,1,1]
%
% Let's break this down. First, filters applied to rows (dimension 1) and
% filters applied to columns (dimension 2) will be handled separately. Filters
% at the top level of the list-of-lists are combined with an AND operation.
% Elements that belong to sub-lists are first combined an OR operation,
% essentially creating a new filter to substitute into the top level of the
% list. This new filter is then combined with AND with other filters of the
% same dimension.
%
% The operation boils down to the following:
%   rowfilter = A & B;
%   G = C | D;
%   colfilter = E & G;
%
% OR operations always precede AND operations. Otherwise, order does not
% matter, since (A & B) == (B & A) and (A | B) == (B | A).
%
% Filters not included in the list of labels will have no effect on the output
% (note that filter F from the example above plays no role).
%
% The function takes no additional arguments, and has only this single mode of
% operation/way of interpretting the syntax of the labels list.
%
% Note that the following will result in an error:
%   composeFilters(FILTERS, {{'A','C'}})
%   Error using composeFilters (line 107)
%   crcox:Cannot OR filters meant for different dimensions.
%
% Chris Cox 22/08/2017
    labels = ascell(labels);

    % metadata.filter points to a structured array of filters.
    % First, force filters to a common orientation.
    for ii = 1:numel(FILTERS)
        FILTERS(ii).filter = forceRowVec(FILTERS(ii).filter);
    end

    % If necessary, combine filters that are meant to be 'OR-ed'
    z = cellfun(@iscell, labels);
    if any(z)
        prototype = [...
            forceRowVec(fieldnames(FILTERS(1)));...
            cell(1,numel(fieldnames(FILTERS(1))))...
        ];
        ix = find(z);
        N = numel(FILTERS);
        FILTERS(end+nnz(z)) = struct(prototype{:});
        to_or = labels(z);
        for i = 1:nnz(z);
            label_set = to_or{i};
            n = numel(label_set);
            if n > 1 % there are multiple filters to OR
                label_set_tagged = cell(1,n+1);
                label_set_tagged(1:n) = label_set;
                label_set_tagged{end} = 'autogen';
                label_comb = strjoin(label_set_tagged, '_');
                filter_selection = selectbyfield(FILTERS(1:N), 'label', label_set);
                if all([filter_selection.dimension]==filter_selection(1).dimension)
                    filter_comb = any(cat(1, filter_selection.filter),1);
                    FILTERS(N+i).label = label_comb;
                    FILTERS(N+i).dimension = filter_selection(1).dimension;
                    FILTERS(N+i).filter = filter_comb;
                    % replace the sublist of filters to OR with the new, combined label
                    % that points to the OR-ed filter.
                    labels{ix(i)} = label_comb;
                else
                    error('crcox:Cannot OR filters meant for different dimensions.')
                end
            else % there aren't multiple filters to OR, we got here by syntactic accident.
                labels(ix(i)) = label_set{1};
            end
        end
    end

    % Then select the filters, and AND them.
    % z = ismember({FILTERS.label}, labels);
    % filters_row = FILTERS(z & [FILTERS.dimension]==1);
    % filters_col = FILTERS(z & [FILTERS.dimension]==2);
    filters_row = selectbyfield(FILTERS, 'label', labels, 'dimension', 1);
    filters_col = selectbyfield(FILTERS, 'label', labels, 'dimension', 2);
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
