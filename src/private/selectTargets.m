function Y = selectTargets(metadata, target, rowfilter)
  Y = {metadata.(target)};
  N = length(Y);
  for i = 1:N
    Y{i} = Y{i}(rowfilter{i});
  end
end
