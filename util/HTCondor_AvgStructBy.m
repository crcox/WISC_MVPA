function [ TA ] = HTCondor_AvgStructBy( results , by, variables, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    p = inputParser();
    addRequired(p, 'results', @isstruct);
    addRequired(p, 'by', @iscellstr);
    addRequired(p, 'variables', @iscellstr);
    addParameter(p, 'concatDim', 1, @isscalar);
    addParameter(p, 'meanDim', [], @isscalar);
    addParameter(p, 'permuteDimsAfterConcat', [], @isvector);
    addParameter(p, 'asTable', false, @islogical);
    parse(p, results, by, variables, varargin{:});
    
    if isempty(p.Results.meanDim)
        meanDim = p.Results.concatDim;
    else
        meanDim = p.Results.meanDim;
    end
    mycat = @(x) catCellsToArray(p.Results.concatDim, p.Results.permuteDimsAfterConcat, x);
    mymean = @(x) mean(x,meanDim);
    T = struct2table(p.Results.results);
    g = findgroups(T(:,p.Results.by));
    TA = unique(T(:,p.Results.by));
    for i = 1:numel(p.Results.variables)
        v = p.Results.variables{i};
        if iscell(T.(v)(1))
            X = splitapply(mycat,T.(v),g);
            x = cellfun(mymean, X, 'UniformOutput', 0);
        else
            x = splitapply(@mean,T.(v),g);
        end
        TA.(v) = x;
    end
    if ~p.Results.asTable
        TA = table2struct(TA);
    end
end

function c = catCellsToArray(dim, perm, x)
    if isempty(perm)
        c = {cat(dim, x{:})};
    else
        c = {permute(cat(dim, x{:}), perm)};
    end
end