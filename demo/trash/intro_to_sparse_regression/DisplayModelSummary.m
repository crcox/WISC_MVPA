function [] = DisplayModelSummary( header, b, bz, yt, ytz, yh, yhz )
    nchar = fprintf('%s\n', header);
    fprintf('%s\n',repmat('=',1,nchar-1));
    fprintf('First 10 weights\n');
    fprintf('----------------\n');
    fprintf('      truth:'); fprintf(' % 8.2f', b(1:10)); fprintf('\n');
    fprintf('  estimated:'); fprintf(' % 8.2f', bz(1:10)); fprintf('\n');
    fprintf('\n');
    fprintf('Accuracy of estimated weights\n');
    fprintf('-----------------------------\n');
    fprintf(' n (true voxels) : % 5d\n', nnz((b~=0) & (bz~=0)));
    fprintf(' n (false voxels): % 5d\n', nnz((b==0) & (bz~=0)));
    fprintf('\n');
    fprintf('Proportion of variance explained\n');
    fprintf('--------------------------------\n');
    fprintf('  R^2 (training set): %.2f\n', corr(yt,ytz)^2);
    fprintf('  R^2 (holdout set) : %.2f\n', corr(yh,yhz)^2);
    fprintf('\n');
end

