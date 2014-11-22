function [RepIndex, groups, group_arr] = define_rep_space(G)
% DEFINE_REP_SPACE Generate all necessary indexes for working with data in a
% replicated space where voxels in multiple groups can be considered as a
% member of each group, and treated like several voxels. Ultimately, the
% weights on several instances of the same voxel in different groups will
% be aggregated.
%
% USAGE:
% [RepIndex, groups, group_arr] = DEFINE_REP_SPACE(G)
%
% INPUTS:
% G:      A Kx1 cell array that contains the column (i.e. voxel) indexes
% in the unreplicated data that belong to each group, where K is the number
% of groups.
%
% OUTPUTS:
% REPINDEX:    A 1xN vector, where N is the number of columns in the
% replicated matrix.  The indexes refer to columns in the unreplicated
% matrix, and when used to select from the unreplicated matrix will produce
% the replicated matrix.  Because the replicated matrix is very large, care
% should be taken to not modify it, otherwise a full copy will be made in
% memory.  If REPINDEX is merely used to index into the unreplicated
% matrix, replication can be achieved by redundant pointers, and not by
% redundent copies.  With large datasets, this is essential to keep the
% problem tractable.
%
% GROUPS:     A 1xN vector, where N has the same meaning as above.  Each
% element in GROUPS is a group label in the replicated space.  The
% replicated matrix has columns ordered by group, so GROUPS will be
% [repmat(1,1,length(G{1})), repmat(2,1,length(G{2})), ... ,
% repmat(K,1,length(G{K}))], where K is the number of groups.
%
% GROUP_ARR: A KxMaxGroupSize matrix, where K means the same as above, and
% MaxGroupSize = max(cellfun(@length,G)).  Rows correspond to groups, and
% elements are indexes into the replicated matrix that would select a
% single group.  If groups are unequal, rows are padded out with N+1, an
% index just outside of the column range of the replicated matrix. This
% dummy index will be handled appropriately at later steps.
%
% DEFINE_REP_SPACE is a memory efficient replacement for makeA_multitask.
%
% Chris Cox | July 24, 2013
    if isrow(G);
        G = G';
    end
    K = length(G);
    MaxGroupSize = uint32(max(cellfun(@length,G)));
    RepSpaceSize = uint32(sum(cellfun(@length,G)));
%% Create RepIndex
    RepIndex = cell2mat(G');

%% Create groups
    temp = G;
    for i = 1:K
        temp{i}(:) = i;
    end
    groups = cell2mat(temp');

%% Create group_arr
    group_arr = nan(K,MaxGroupSize);
    a = 1;
    for i = 1:K
        glen = length(G{i});
        b = a+glen-1;
        group_arr(i,1:glen) = a:b;
        a = b+1;
    end
end
