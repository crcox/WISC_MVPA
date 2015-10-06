function C = ascell(X)
  if ~iscell(X)
    C = {X};
  else
    C = X;
  end
end
