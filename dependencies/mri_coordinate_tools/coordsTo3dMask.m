function [mask,dcoords] = coordsTo3dMask(coords)
  dcoords = voxelspacing(coords);
  mcoords = min(coords);
  ijk = roundto(bsxfun(@minus, coords, mcoords),dcoords);
  ijk = roundto(bsxfun(@rdivide, ijk, dcoords),1)+1;
  ind = sub2ind(max(ijk), ijk(:,1), ijk(:,2), ijk(:,3));
  mask = false(max(ijk));
  mask(ind) = true;
end
