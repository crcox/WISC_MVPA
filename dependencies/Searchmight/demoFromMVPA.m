function [] = demoFromMVPA()

% adapted from  
%
%    TUTORIAL_EASY.HTM 
%    This is the sample script for the Haxby et al. (Science, 2001) 8-
%    categories study. See the accompanying TUTORIAL_EASY.HTM, the
%    MVPA manual (MANUAL.HTM) and then TUTORIAL_HARD.HTM.
%
% requires SearchmightToolbox in your path, so
% - addpath(<path to SearchmightToolbox)
% - setupPathsSearchmightToolbox
%
% and also the path you'd need in order to run the MVPA toolbox tutorial
%

% start by creating an empty subj structure
subj = init_subj('haxby8','tutorial_subj');

%%% create the mask that will be used when loading in the data
subj = load_afni_mask(subj,'VT_category-selective','mask_cat_select_vt+orig');

% now, read and set up the actual data. load_AFNI_pattern reads in the
% EPI data from a BRIK file, keeping only the voxels active in the
% mask (see above)
for i=1:10
  raw_filenames{i} = sprintf('haxby8_r%i+orig',i);
end
subj = load_afni_pattern(subj,'epi','VT_category-selective',raw_filenames);

% initialize the regressors object in the subj structure, load in the
% contents from a file, set the contents into the object and add a
% cell array of condnames to the object for future reference
subj = init_object(subj,'regressors','conds');
load('tutorial_regs');
subj = set_mat(subj,'regressors','conds',regs);
condnames = {'face','house','cat','bottle','scissors','shoe','chair','scramble'};
subj = set_objfield(subj,'regressors','conds','condnames',condnames);

% store the names of the regressor conditions
% initialize the selectors object, then read in the contents
% for it from a file, and set them into the object
subj = init_object(subj,'selector','runs');
load('tutorial_runs');
subj = set_mat(subj,'selector','runs',runs);

%  condnames          1x8                   566  cell                
%  raw_filenames      1x10                  882  cell                
%  regs               8x1210              77440  double              
%  runs               1x1210               9680  double              
%  subj               1x1               7237645  struct    

% get data into a #TRs x #voxels matrix, and also mask

TRdata = subj.patterns{1}.mat';
[nTRs,nVoxels] = size(TRdata);

mask = subj.masks{1}.mat;
[dimx,dimy,dimz] = size(mask);

% create labels for each TR (column vectors)

TRlabels      = sum(regs .* repmat((1:8)',1,nTRs),1)';
TRlabelsGroup = runs';

%
% turn each block into an example, and convert labels
%

% binary mask for beginning and end of blocks (task or fixation)
TRmaskBlockBegins = ([1;diff(TRlabels)] ~= 0); 
TRmaskBlockEnds   = ([diff(TRlabels);1] ~= 0);

% average blocks (we will get rid of 0-blocks afterwards)
% (silly average of all images, without thinking of haemodynamic response)
% to convert them into examples (and create labels and group labels)

% figure out how many blocks and what TRs they begin and end at
nBlocks = sum(TRmaskBlockBegins);
blockBegins = find(TRmaskBlockBegins);
blockEnds   = find(TRmaskBlockEnds);

% create one example per block and corresponding labels
labels      = zeros(nBlocks,1); % condition
labelsGroup = zeros(nBlocks,1); % group (usually run)
examples    = zeros(nBlocks,nVoxels); % per-block examples

for ib = 1:nBlocks
  range = blockBegins(ib):blockEnds(ib);
  examples(ib,:)  = mean(TRdata(range,:),1);
  labels(ib)      = TRlabels(blockBegins(ib));
  labelsGroup(ib) = TRlabelsGroup(blockBegins(ib));
end

% nuke examples with label 0
indicesToNuke = find(labels == 0);
examples(indicesToNuke,:)  = [];
labels(indicesToNuke)      = [];
labelsGroup(indicesToNuke) = [];

% create a meta structure (see README.datapreparation.txt and demo.m for more details about this)
meta = createMetaFromMask(mask);

%
% run a fast classifier (see demo.m for more details about computeInformationMap)
%

classifier = 'gnb_searchmight'; % fast GNB

[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'searchlight', ...
                                meta.voxelsToNeighbours,meta.numberOfNeighbours);

%
% quick plot of the results
%

clf; nrows = ceil(sqrt(dimz)); ncols = nrows;

volume = repmat(NaN,[dimx dimy dimz]);
% place accuracy map in a 3D volume, using the vectorized indices of the mask in meta
volume(meta.indicesIn3D) = am;

for iz = 1:dimz
  subplot(nrows,ncols,iz);
  imagesc(volume(:,:,iz)',[0 0.5]); axis square;
  set(gca,'XTick',[]); set(gca,'YTick',[]);
  if iz == 1; hc=title('accuracy map for mask'); set(hc,'FontSize',8); end
  if iz == dimz; hc=colorbar('vert'); set(hc,'FontSize',8); end
end

