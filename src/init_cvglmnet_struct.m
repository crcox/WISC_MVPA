function y = init_cvglmnet_struct(ncv, varargin)
    if numel(varargin) == 0
        dims = 1;
    else
        if prod(varargin{1}) > max(varargin{1})
            dims = varargin{1};
        else
            dims = cell2mat(varargin);
        end
    end
    x = cvglmnet(magic(100),repmat([0;1],50,1),'binomial');
    f = fieldnames(x);
    X = cell2struct(repmat({[]}, numel(f), 1), f);
    f = fieldnames(x.glmnet_fit);
    X.glmnet_fit = repmat(cell2struct(repmat({[]}, numel(f), 1), f), 1,ncv);
    y = repmat(X, dims);
end