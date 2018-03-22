function [ TopN_by_group, ix ] = hyperband_pick_top_n( ModelInstances, n, maximize)
%HYPERBAND_PICK_TOP_N Select models that survive current hyperband round
%   Provided a set of Model instances (bundled in ModelContainers)
    if nargin < 3
        maximize = false;
    end
    if maximize
        direction = 'descend';
    else
        direction = 'ascend';
    end
    
    
    regularization = ModelInstances(1).regularization;
    switch upper(regularization)
        case 'SOSLASSO'
            testError = cellfun(@(x) mean(x.Model.testError), num2cell(ModelInstances))';
            t0 = struct2table(rmfield(ModelInstances,{'Model','cvholdout','sim_source','sim_metric','data','subject','G'})');
        case {'LASSO','RIDGE'}
            testError = cellfun(@(x) x.Model.testError, num2cell(ModelInstances))';
            t0 = struct2table(rmfield(ModelInstances,{'Model','cvholdout','sim_source','sim_metric','G'})');
        otherwise
            testError = cellfun(@(x) x.Model.testError, num2cell(ModelInstances))';
            t0 = struct2table(rmfield(ModelInstances,{'Model','cvholdout'})');     
    end
    [t1,~,g0] = unique(t0);
    g0_max = max(g0);
    t1.cvmap = cell(size(t1,1),1);
    for i = 1:g0_max
        t1.cvmap{i} = find(g0 == i);
    end
    t1.testError = grpstats(testError, g0);
    t1.rank = zeros(size(t1,1),1);
    t2 = t1;
    switch upper(regularization)
        case {'SOSLASSO','LASSO'}
            try
                t2.lambda = [];
            catch
                t2.lamL1 = [];
            end
            try
                t2.alpha = [];
            catch
                t2.lamSOS = [];
            end
        case {'RIDGE'}
            t2.lamL2 = [];
            
        otherwise
            t2.lambda = [];
            t2.lambda1 = [];
    end

    t2.configID = [];
    t2.testError = [];
    t2.cvmap = [];
    [~,~,g1] = unique(t2);
    g1_max = max(g1);
    for i = 1:g1_max
        tmp = t1(g1==i,:);
        [~,ix] = sort(tmp.testError, direction);
        t1(g1==i,:) = tmp(ix,:);
        t1.rank(g1==i) = 1:numel(ix);
    end
    ix = cell2mat(t1.cvmap(t1.rank<=n));
    TopN_by_group = ModelInstances(ix);
end

