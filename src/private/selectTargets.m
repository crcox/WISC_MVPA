function Y = selectTargets(metadata, target, rowfilter)
  Y = cell(1, numel(metadata));
  if numel(metadata) == 1;
    rowfilter = {rowfilter};
  end
  for i = 1:numel(metadata);
    M = metadata(i);
    z = strcmp(target, {M.target.label});
    Y{i} = M.targets(z).target(rowfilter{i});
  end
end
