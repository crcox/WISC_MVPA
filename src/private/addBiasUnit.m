function X = addBiasUnit(X)
  X = ascell(X);
  N = length(X);
  for i = 1:N
    X{i} = [X{i}, ones(size(X{i},1),1)];
  end
  X = uncell(X);
end
