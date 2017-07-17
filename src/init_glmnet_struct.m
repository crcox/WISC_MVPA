function y = init_glmnet_struct(varargin)
    if nargin == 0
        dims = 1;
    else
        if prod(varargin{1}) > max(varargin{1})
            dims = varargin{1};
        else
            dims = cell2mat(varargin);
        end
    end
    x = glmnet(magic(10),repmat([0;1],5,1),'binomial');
    f = fieldnames(x);
    x = cell2struct(repmat({[]}, numel(f), 1), f);
    y = repmat(x, dims);
end