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
    hyperparams_c = cellfun(@(x) strjoin({x,'c'},'_'),hyperparams,'UniformOutput',false);
    for i = 1:numel(hyperparams)
        try
            T.(hyperparams_c{i}) = categorical(T.(hyperparams{i}));
        catch
            T.(hyperparams_c{i}) = categorical(cell2mat(T.(hyperparams{i})));
        end
    end
    
    for i = 1:numel(p.Results.by)
        try
            T.(p.Results.by{i}) = categorical(T.(p.Results.by{i}));
        catch
            T.(p.Results.by{i}) = categorical(cell2mat(T.(p.Results.by{i})));
        end
    end
    
    [G, Tavg, ia] = AssignGroupLabels(T(:,[hyperparams_c,p.Results.by]));

    for i = 1:numel(hyperparams)
        Tavg.(hyperparams{i}) = T.(hyperparams{i})(ia);
    end
    
    varsToAverage = [{objective},p.Results.extras];
    for i = 1:numel(varsToAverage)
        f = varsToAverage{i};
        Tavg.(f) = ApplyByGroup(@mean, T.(f), G);
    end

    G = AssignGroupLabels(Tavg(:,p.Results.by));
    
    varsToReport = [{objective},p.Results.extras,hyperparams,p.Results.by];
    Tbest = OptByGroup(objective, Tavg(:,varsToReport), G, 'minimize', p.Results.minimize);
    
% Solution A
% ==========
%     Tbest = array2table(X,'VariableNames',varsToReport);
%     Tbest = cat(1, X{:});

% Solution B
% ==========
%     for i = 1:numel(varsToReport)
%         f = varsToReport{i};
%         if isnumeric(T.(f))
%             Tbest.(f) = cell2mat(X(:,i));
%         else
%             Tbest.(f) = X(:,i);
%         end
%     end
end

% 
% function [X,I] = whichmin( varargin )
%     objective = varargin{1};
%     X = cell(size(varargin));
%     [~,I] = min(objective);
%     for i = 1:numel(varargin)
%         X{i} = varargin{i}(I,:);
%     end
% end

function y = OptByGroup(objective, x, g, varargin)
    p = inputParser();
    addRequired(p, 'objective');
    addRequired(p, 'x');
    addRequired(p, 'g');
    addParameter(p, 'minimize', true);
    parse(p, objective, x, g, varargin{:});

    if iscell(x)
        x = cell2mat(x);
    end
    
    if isnumeric(objective)
        z = (1:size(x,2)) == objective;
    end
        
    if istable(x)
        tin = true;
        labs = x.Properties.VariableNames;
        q = x.(objective);
    else
        tin = false;
        q = x(:, z);
    end

    if p.Results.minimize
        ind = accumarray(g, q, [], @whichmin);
    else
        ind = accumarray(g, q, [], @whichmax);
    end
    
%     y = zeros(max(g), size(x,2));
    y = x(1:max(g),:);
    for i = 1:max(g)
        tmp = x(g == i, :);
        y(i,:) = tmp(ind(i), :);
    end
    
    function ind = whichmax( x )
        [~,ind] = max(x);
    end

    function ind = whichmin( x )
        [~,ind] = min(x);
    end 
end

function y = ApplyByGroup(func, x, g)
    if iscell(x)
        x = cell2mat(x);
    end
        
    if exist('splitapply','builtin') == 5
        y = splitapply(func, x, g);
    else
        if istable(x)
            tin = true;
            labs = x.Properties.VariableNames;
            x = table2array(x);
        else
            tin = false;
        end
        
        y = accumarray(g, x, [], func);
        
        if tin
            y = array2table(y, 'VariableNames', labs);
        end
    end
end

function [g,avg,ia] = AssignGroupLabels(x)
    if exist('findgroups', 'builtin') == 5
        g = findgroups(x);
        [avg,ia] = unique(x);
    else
        [avg,ia,g] = unique(x);
    end
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
