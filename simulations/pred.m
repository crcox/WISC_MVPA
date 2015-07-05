function Yz = pred(X,bz)
  Yz = zeros(size(X{1},1), numel(X));
  for i = 1:numel(X)
    Yz(:,1) = X{i} * bz(:,i);
  end
end