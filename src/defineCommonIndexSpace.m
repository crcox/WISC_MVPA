function [I,J,K] = defineCommonIndexSpace(XYZ,SharedSpaceVoxelSize,GroupSize,GroupShift)
%% Find the range of each dimension that contains all coordinates.
    minXYZ = min(XYZ);
    maxXYZ = max(XYZ);
    maxIJK = (maxXYZ - minXYZ)+ 1;
    maxIJK = bsxfun(@rdivide,maxIJK,SharedSpaceVoxelSize);
    I = maxIJK(1);
    J = maxIJK(2);
    K = maxIJK(3);
    clear minXYZ maxXYZ maxIJK;

    %% Inflate slightly, as needed, to accomodate group size (and GroupShift).
    I = I + (GroupShift(1)-mod(I - GroupSize(1),GroupShift(1)));
    J = J + (GroupShift(2)-mod(J - GroupSize(2),GroupShift(2)));
    K = K + (GroupShift(3)-mod(K - GroupSize(3),GroupShift(3)));
end
