function y = nrsa_corr(x,p)
% NRSA_CORR Correlation-base performance metric for Network RSA
%
% Computes the Spearman correlation between a target similarity matrix and
% a predicted similarity matrix (lower triangles).
%
% INPUT
% x : the target matrix
% p : the predicted matrix
%
% OUTPUT
% y : error term (the loss value)
    z = tril(true(size(x)),-1);
    y = corr(x(z),p(z),'type','Spearman');
end