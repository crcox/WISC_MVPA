function M = uncell(C)
  % If C is a cell with one element, extract that element and descard the cell
  % casing.
  M = C;
  if iscell(C)
    if length(C) == 1
      M = C{1};
    end
  end
end
