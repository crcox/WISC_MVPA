function X = proxSortedL1L2(Y,lambda)


y = sqrt(sum(Y.^2,2));
r = size(Y,2);

x = proxSortedL1(y,lambda);

X = Y .* repmat(x./(y + realmin),[1,r]);

