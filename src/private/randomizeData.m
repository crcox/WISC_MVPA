function X = randomizeData(X) %#ok<DEFNU>
  X = ascell(X);
  N = length(X);
  for i = 1:N
    X{i} = randn(size(X{i}));
  end
  if N == 1
    X = X{i};
  end
  X = uncell(X);
end
