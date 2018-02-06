function S = cosinesimilarity(X);
  [n,d] = size(X);
  X = X./(sqrt(sum(X.^2,2))*ones(1,d));
  S = X*X';
end
