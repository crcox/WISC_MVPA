function [X,GroupInfo] = prep_data(X,metadata,params)
	%% Parameters
	nSubjects  = length(X);
	GroupSize  = params.GroupSize;
	GroupShift = params.GroupShift;
	XYZ_tlrc   = cell(nSubjects,1);

	%% Round coordinates to nearest multiple of SharedSpaceVoxelSize.
	for ss = 1:nSubjects
		XYZ_tlrc{ss} = roundto(metadata(ss).xyz_tlrc, params.SharedSpaceVoxelSize);
	end

	%% Concatenate all coordinates
	XYZ_tlrc_mat = cell2mat(XYZ_tlrc);

	%% Find the range of each dimension that contains all coordinates.
	[I,J,K] = defineCommonIndexSpace(XYZ_tlrc_mat, ...
			params.SharedSpaceVoxelSize, ...
			params.GroupSize, ...
			params.GroupShift);

	GroupInfo.ijk_tlrc_range = [I J K];

	%% MAPPING FROM ijk ---> n = i + I(j-1) + IJ(k-1)
	IND_tlrc = cell(nSubjects,1);
	IJK_tlrc = cell(nSubjects,1);
	for s = 1:nSubjects
		[IJK_tlrc{s}, IND_tlrc{s}] = xyz2ijk(XYZ_tlrc{s}, [I,J,K], params.SharedSpaceVoxelSize);
	end
	IND_tlrc_vec = cell2mat(IND_tlrc);
	% Breakdown:
	% 1. Transformation to 1 based IJK indexes in TLRC space: Take
	% coordinates for all subjects as one large matrix, subtract the
	% minimums, and add 1.
	% 2. Transform those IJK indexes to literal column indexes. Retain only the
	% unique ones, and store as an array of unsigned 32 bit integers.

	% Construct a mask for the whole common space, to identify only those
	% points the contain data.
	NZ_tlrc = false(1,I*J*K);
	NZ_tlrc(IND_tlrc_vec) = true;
	GroupInfo.NZ_tlrc = NZ_tlrc;

	% MAKE NEW DATA MATRICES
	fprintf('\n');
	SubjectIndexes = cell(nSubjects,1);
	SubjectMask = false(sum(NZ_tlrc),nSubjects);
	% Subject Mask and Subject Indexes work together. The Subject Mask makes
	% it clear which weights/voxels belong to each subject. The Subject Index
	% will put them back in that subject's original order.
	for ss = 1:nSubjects
		fprintf('% 2d',ss)
		[~,ind_tlrc] = xyz2ijk(XYZ_tlrc{ss}, [I,J,K], params.SharedSpaceVoxelSize);
		[~,ix] = sort(ind_tlrc);
		[~,SubjectIndexes{ss}] = sort(ix);
		% At this line, the dataset is huge:
		funcData_tlrc = zeros(size(X{ss},1),I*J*K);
		funcData_tlrc(:,ind_tlrc) = X{ss};
		funcData_tlrc = sparse(funcData_tlrc(:,NZ_tlrc));
		SubjectMask(:,ss) = any(funcData_tlrc);
		X{ss} = funcData_tlrc;
		clear ijk xyz_tlrc;
	end
	fprintf('\n')

	%% Define Range
	if all(XYZ_tlrc_mat(:,1)==1)
		irange = 1;
	else
		irange = 1:GroupShift(1):(I-GroupSize(1));
	end
	if all(XYZ_tlrc_mat(:,2)==1)
		jrange = 1;
	else
		jrange = 1:GroupShift(2):(J-GroupSize(2));
	end
	if all(XYZ_tlrc_mat(:,3)==1)
		krange = 1;
	else
		krange = 1:GroupShift(3):(K-GroupSize(3));
	end

	%% Generate Grid of Corners
	[KG,JG,IG] = ndgrid(krange,jrange,irange);
	IJK_corners = [IG(:), JG(:), KG(:)]; % Grid of group corners.

	%% Generate Grid-patch to make groups
	[IG,JG,KG]= ndgrid(0:(GroupSize(1)-1),0:(GroupSize(2)-1),0:(GroupSize(3)-1));
	GroupExtent = [IG(:), JG(:), KG(:)]; % Grid of distances from corner.

	%% Make G (full)
	% At this point, the indexes will refer into the FULL common space.
	% They will need to be adjusted later so that they actually point to
	% columns in X{i}, which have empty voxels removed.
	M = length(IJK_corners);
	G = cell(M,1);
	for i = 1:(M)
		G_ijk = bsxfun(@plus, IJK_corners(i,:), GroupExtent);
		temp = sub2ind([I,J,K], G_ijk(:,1),G_ijk(:,2),G_ijk(:,3))';
		G{i} = uint32(temp(NZ_tlrc(temp)));
	end

	% This will just stick any voxels that aren't yet in a group into a group
	% together.
	CatchallGroup = setdiff(IND_tlrc_vec,cell2mat(G'));
	if ~isempty(CatchallGroup)
		G{M+1} = uint32(CatchallGroup(NZ_tlrc(CatchallGroup)));
	end
	clear IND_tlrc_vec;
	%% Remove Empty Groups
	G = G(~cellfun('isempty',G));

	%% Map to small space
	GroupInfo.G_oldspace = G;
	M = length(G);
	G = cell(M,1);

	nth_NZ_TLRC = repmat(uint32(0),1,length(NZ_tlrc));
	nth_NZ_TLRC(NZ_tlrc) = 1:sum(NZ_tlrc);
	for ii = 1:M
			G{ii} = nth_NZ_TLRC(GroupInfo.G_oldspace{ii});
	end
	GroupInfo.G = G;
	clear G;

	%% Define the replicated, non-overlapping space.
	[GroupInfo.RepIndex, GroupInfo.groups, GroupInfo.group_arr] = define_rep_space(GroupInfo.G);
	GroupInfo.SubjectIndexes = SubjectIndexes;
	GroupInfo.SubjectMask = SubjectMask;
end
