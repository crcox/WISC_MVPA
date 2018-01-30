function [ TopN_by_group, ix ] = hyperband_pick_top_n( AdlasInstances, n, maximize)
%HYPERBAND_PICK_TOP_N Select models that survive current hyperband round
%   Provided a set of Adlas instances (bundled in AdlasContainers)
    if nargin < 3
        maximize = false;
    end
    if maximize
        direction = 'descend';
    else
        direction = 'ascend';
    end

    testError = cellfun(@(x) x.Adlas.testError, num2cell(AdlasInstances))';
    t0 = struct2table(rmfield(AdlasInstances,{'Adlas','cvholdout'})');
    [t1,~,g0] = unique(t0);
    g0_max = max(g0);
    t1.cvmap = cell(size(t1,1),1);
    for i = 1:g0_max
        t1.cvmap{i} = find(g0 == i);
    end
    t1.testError = grpstats(testError, g0);
    t1.rank = zeros(size(t1,1),1);
    t2 = t1;
    t2.lambda = [];
    t2.lambda1 = [];
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
    TopN_by_group = AdlasInstances(ix);
end

