function C = selectCoordinates(metadata, orientation, type, colfilter)
  C = cell(1, numel(metadata));
  for i = 1:numel(metadata);
    c = metadata(i).coords;
    z = all([strcmpi(orientation,{c.orientation});strcmpi(table,{c.type})]);
    C{i} = metadata(i).coords(z).(type)(colfilter{i});
  end
end
