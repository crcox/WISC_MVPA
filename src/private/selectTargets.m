function Y = selectTargets(metadata, target, rowfilter)
  Y = cell(1, numel(metadata));
  for i = 1:numel(metadata);
    M = metadata(i);
    z = strcmp(target, {M.targets.label});
    Y{i} = M.targets(z).targets(rowfilter{i});
  end
end
