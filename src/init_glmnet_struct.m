function y = init_glmnet_struct(opts)
    if nargin == 0
        opts = glmnetSet();
    else
        opts = glmnetSet(opts);
    end
    x = glmnet(magic(10),repmat([0;1],5,1),'binomial',opts);
    f = fieldnames(x);
    x = cell2struct(repmat({[]}, numel(f), 1), f);
    y = repmat(x, 1); % fix later maybe?
end