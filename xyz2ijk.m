function [ijk, ind] = xyz2ijk(xyz, commonspacedims, voxelsize)
    ijk = bsxfun(@minus,xyz,min(xyz));
    ijk = bsxfun(@rdivide, ijk, voxelsize) + 1;
    if any(size(ijk) < 2)
        ind = ijk;
    else
        ijksplat = mat2cell(ijk,size(ijk,1), ones(1,size(ijk,2)));
        ind = uint32(sub2ind(commonspacedims, ijksplat{:}));
    end
end