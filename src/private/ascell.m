function C = ascell(X)
  if ~iscell(X) && ~isstring(X)
    C = {X};
  else
    C = X;
  end
end
