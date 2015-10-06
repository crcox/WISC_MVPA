function X = shuffleData(X,rowfilter) %#ok<DEFNU>
  X = ascell(X);
  N = length(X);
  n = length(rowfilter{1});
  shuffledIndex = randperm(n);
  for i = 1:N
    ix   = rankind(shuffledIndex(rowfilter{i}));
    X{i} = X{i}(ix,:);
  end
  X = uncell(X);
end
