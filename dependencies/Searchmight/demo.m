%
% this function creates a synthetic dataset and shows various examples of
% how to call the code (scroll down)
%

function [] = demo()

fprintf('Please do not run this directly, just cut and paste code from the commented code below or copy this function into a new one and edit that instead. Thank you!\n'); return;

%
% this is a mock dataset, a cubic brain floating in a 11x11x12 volume (a blockhead?)
%

% create a 3D binary mask 

dimx = 11; dimy = 11; dimz = 12;
mask = zeros(dimx,dimy,dimz);
xrange = 2:10; yrange = 2:10;
for iz = 2:(dimz-1); mask(xrange,yrange,iz) = 1; end

% create the "meta" neighbourhood structure

meta = createMetaFromMask(mask);
nVoxels = length(meta.indicesIn3D);

%% generate data

% create per-class activity templates, one activation spot per class
nClasses = 2;

voxelMeans3D = zeros(dimx,dimy,dimz,nClasses);
voxelMeans3D(5,5,3,1) = 1;
voxelMeans3D(5,5,5,1) = 1;
voxelMeans3D(5,5,7,2) = 1;
voxelMeans3D(5,5,9,2) = 1;

for ic = 1:nClasses
  tmp = voxelMeans3D(:,:,:,ic);
  voxelMeans{ic} = tmp(meta.indicesIn3D)';
end

% create examples (4 groups, think of these as runs)

nGroups = 4;
nPerClassPerGroup = 25; npcpg = nPerClassPerGroup;

n = nGroups * nPerClassPerGroup * 2;
m = nVoxels;

examples    = zeros(n,m); % patterns
labels      = zeros(n,1); % labels of each image
labelsGroup = zeros(n,1); % cross-validation labels ("runs", say)

idx = 1;
for ig = 1:nGroups
  for ic = 1:nClasses
    range = idx:(idx+npcpg-1);
    
    examples(range,:) = randn(npcpg,m) * 0.5;
    examples(range,:) = examples(range,:) + repmat(voxelMeans{ic},npcpg,1);
    
    labels(range)      = ic;
    labelsGroup(range) = ig;
    idx = idx + npcpg;
  end
end

labels = labels + 1;

%
% if you want to see what an example looks like
%

example = examples(1,:);

% place the example in a 3D volume
volume = repmat(NaN,meta.dimensions); 
volume(meta.indicesIn3D) = example;

% plot it in 3 x 4
nr = 3; nc = 4; clf;
for iz=1:12; subplot(nr,nc,iz); imagesc(volume(:,:,iz),[-1 1]); axis square; end


%
% a few example command lines with the various types of classifier
% (am is an accuracy map, pm is the corresponding analytical p-value map,
%  plot code for results at the bottom)
%

% single voxel
classifier = 'gnb_pooled';
[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier);

% searchlight, gaussian bayesian classifiers
classifier = 'gnb_pooled';
classifier = 'lda_shrinkage';
classifier = 'gnb_searchmight'; % fast GNB

[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'searchlight', ...
                                meta.voxelsToNeighbours,meta.numberOfNeighbours);


%
% what can you do with these?
%

%% quick plot of any accuracy map

% place the accuracy map in a 3D volume using the meta structure
volume = repmat(NaN,meta.dimensions); 
volume(meta.indicesIn3D) = am;

% plot proper
nr = 3; nc = 4; clf;
for iz=1:12; subplot(nr,nc,iz); imagesc(volume(:,:,iz),[0 1]); end

% threshold a p-value map using a clone of Tom Nichols' FDR function FDR.m:
q = 0.01;
[thresholdID,thresholdN] = computeFDR(pm,q);  % produces two thresholds, use pN;
binaryMap = pm <= thresholdN; % thresholded map

%
% the previous p-values were analytical, you can also get permutation p-values
%

classifier = 'gnb_searchmight'; % fast GNB, other classifiers will take a lot longer...

[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'searchlight', ...
                                meta.voxelsToNeighbours,meta.numberOfNeighbours,'testToUse','accuracyOneSided_permutation',1000);

%
% more classifiers
%

% searchlight, nearest-something classifiers

classifier = 'knn'; classifierParameters={'k',0}; % nearest class mean
classifier = 'knn'; classifierParameters={'k',1}; % nearest neighbour
classifier = 'knn'; classifierParameters={'k',3}; % 3-nearest-neighbours

% all of these can work with different similarity/distance measures
similarityMeasure = 'correlation';
similarityMeasure = 'euclidean';
classifier = 'knn'; classifierParameters={'k',1,'similarityMeasure',similarityMeasure};

[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'classifierParameters',classifierParameters,'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours);

% searchlight, SVM classifiers

% all with LIBSVM defaults (C=1, gamma=1/#features
classifier = 'svm_linear';
classifier = 'svm_rbf';
classifier = 'svm_quadratic';

[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'searchlight', ...
                                meta.voxelsToNeighbours,meta.numberOfNeighbours);

% specify parameters
classifier = 'svm_linear'; classifierParameters = {'lambda',1};
classifier = 'svm_rbf'; classifierParameters = {'lambda',1,'gamma',1};

[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'classifierParameters',classifierParameters,'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours);

% use cross-validation within the training set to pick parameters (time consuming)
classifier = 'svm_linear'; classifierParameters = {'lambda','crossValidation'};
classifier = 'svm_rbf'; classifierParameters = {'lambda','crossValidation','gamma','crossValidation'};

[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'classifierParameters',classifierParameters,'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours);

%
% extract additional information per searchlight (comes as a field in extraReturns)
%

classifier = 'lda_shrinkage';

% covariance matrices (supported in LDA or QDA only)
[am,pm,extraReturns] = computeInformationMap(examples,labels,labelsGroup,classifier, ...
                                             'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours,'storeSearchlightCovarianceMatrices');


% they are stored as the field
% extraReturns.searchlightCovarianceMatrices: [27x27x810 double]
% (one 27x27 matrix per radius-1 searchlight for each of 810 voxels)

% confusion matrices (supported in all classifiers)
[am,pm,extraReturns] = computeInformationMap(examples,labels,labelsGroup,classifier, ...
                                             'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours,'storeLocalConfusionMatrices');

% they are stored as the field
% extraReturns.localConfusionMatrices: [2x2x810 double]
% (one #classes x #classes matrix per radius-1 searchlight for each of 810 voxels)

%
% other things
%

% compute autocorrelation between successive examples at each voxel
% (as suggested in the paper as a coarse test of whether examples from consecutive trials/blocks
% are independent (lag = 1 for consecutive examples, lag = 2 for examples with one in-between)

[autocorrelation] = computeAutocorrelation(examples,[1 2]);

clf;
subplot(1,2,1);
hist(autocorrelation(1,:),20);title('autocorrelation at lag 1');axis([-1 1 0 Inf]);
subplot(1,2,2);
hist(autocorrelation(2,:),20);title('autocorrelation at lag 2');axis([-1 1 0 Inf]);

