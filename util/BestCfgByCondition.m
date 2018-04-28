function [ Tbest, Tavg ] = BestCfgByCondition( tune_table, hyperparams, objective, varargin)
%BESTCFGBYCONDITION Reports a table of configurations based on parameter
%search.
%   Inputs:
%   <required>
%   tune_dir     A directory filed with a numbered subdirectory for each
%                job the the hyperparameter search was split into. Each
%                subdirectory should contain, at least, a results.mat (as
%                output from WholeBrain_MVPA) and a params.json file (input
%                cfg for WholeBrain_MVPA).
%   hyperparams  A cellstr contains a list of parameters that were being
%                tuned. For example, lamSOS, lamL1, and diameter for SOS
%                Lasso.
%   objective    A char (string) indicating which variable should be
%                minimized (or maximized, if minimize == false) when
%                selecting the best cfg.
%   <optional>
%   extras       A cellstr containing a list of variables that should be
%                reported along with the best cfg and the objective value.
%   by           A cellstr containing one or more variables to group by
%                when picking best configs. Useful if multiple conditions
%                were tuned at the same time. (To use without extras, set
%                extras to {}).
%   <keyword args>
%   minimize     If true (default), the objective value will be minimized
%                when reporting the best configs. If false, it will be
%                maximized.
%
    p = inputParser();
    addRequired(p, 'tune_table', @istable);
    addRequired(p, 'hyperparams', @iscellstr);
    addRequired(p, 'objective', @ischar);
    addOptional(p, 'extras', {}, @iscellstr);
    addOptional(p, 'by', {}, @iscellstr);
    addParameter(p, 'minimize', true, @(x) islogical(x) || any(x==[0,1]));
    parse(p, tune_table, hyperparams, objective, varargin{:});
    
    T = p.Results.tune_table;
%     T = handle_empty_rows(T,hyperparams,p.Results.by);      
    
    try
        G = findgroups(T(:,[hyperparams,p.Results.by]));
        Tavg = unique(T(:,[hyperparams,p.Results.by]));
    catch ME
        [Tavg,~,G] = unique(T(:,[hyperparams,p.Results.by]));
    end
    varsToAverage = [{objective},p.Results.extras];
    for i = 1:numel(varsToAverage)
        f = varsToAverage{i};
        Tavg.(f) = splitapply(@mean, T.(f), G);
    end

    try
        G = findgroups(Tavg(:,p.Results.by));
        Tbest = unique(Tavg(:,p.Results.by));
    catch ME % older versions do not have find groups
        [Tbest,~,G] = unique(Tavg(:,p.Results.by));
    end
    
    varsToReport = [{objective},p.Results.extras,hyperparams];
    if p.Results.minimize
        whichmin_obj = @(x) whichmin(x);
        X = splitapply(whichmin_obj, Tavg{:,varsToReport}, G);
    else
        whichmax_obj = @(x) whichmax(objective, x);
        X = splitapply(whichmax_obj, Tavg(:,varsToReport), G);
    end
%     Tbest = array2table(X,'VariableNames',varsToReport);
%     Tbest = cat(1, X{:});
    for i = 1:numel(varsToReport)
        f = varsToReport{i};
        Tbest.(f) = X(:,i);
    end
end

function [X,I] = whichmin( x )
    [~,I] = min(x(:,1));
    X = x(I,:);
end

function [X,I] = whichmax( objective, x )
    [~,I] = max(x.(objective));
    X = x(I,:);
end

% function T = handle_empty_rows(T,hyperparams,by)
%     z = cellfun(@isempty, T.(hyperparams{1}));
%     if any(z)
%         T(z,:) = [];
%         f = [hyperparams,by];
%         for i = 1:numel(f)
%             x = T.(f{i})(1);
%             if iscell(x) && isnumeric(x{1}) && isscalar(x{1})
%                 try
%                     T.(f{i}) = cell2mat(T.(f{i}));
%                 catch ME
%                     disp('Expected a scalar, but did not receive one.');
%                     rethrow(ME);
%                 end
%             end
%         end
%     end
% end
