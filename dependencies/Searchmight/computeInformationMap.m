%
% Compute an information map based on classification accuracy at each voxel or a region around
% it, using cross-validation.
%
% In:
% - examples (#examples x #features array)
% - labels   (#examples x 1 vector)
% - group labels (#examples x 1 vector) (could be run labels, for instance, used for cross-validation)
% - classifier to use
%   - gaussian
%     - 'gnb' - a gaussian naive bayes (assuming same variance for both classes)
%     - 'gnb_searchmight' - same, MEX version that runs fast
%     - 'lda' - pure LDA
%     - 'lda_ridge' - LDA with a small ridge term
%     - 'lda_shrinkage' - LDA with a shrinkage estimator for the covariance matrix (Strimmer
%     2005)
%     - 'qda_shrinkage' - QDA with same, given that lda_shrinkage is better than the others might
%     as well use it
%   - SVM 
%     - 'svm_linear'
%     - 'svm_quadratic'
%     - 'svm_sigmoid'
%     - 'svm_rbf'
%
% - test to use (see Pereira/Mitchell/Botvinick 2009)
%     - 'accuracyOneSided_analytical' - a simple one-sided binomial test
%     - 'accuracyOneSided_permutation',<nPermutations> - a permutation test
%       (with gnb_searchmight it's feasible to do 100K permutations overnight, all others are way slower)
%
% - 'searchlight',voxelsToNeighbours,numberOfNeighbours
%
%   in order to use any searchlight classifier you need extra information about which voxels are
%   adjacent to which other voxels in space. This code uses two matrices that encode the
%   relationships, instead of picking one specific neighbourhood type (such as a cube or a
%   sphere). Look at createMetaFromMask.m to see how to take a binary 3D matrix containing a
%   mask of which voxels are in the brain and produce a structure which has
%   voxelsToNeighbours,numberOfNeighbours as fields (or README.datapreparation.txt)
%
% Out:
% - a map of accuracy at each voxel or region around it
% - a map of p-values of the test of the null hypothesis at each voxel 
%
% Examples
% 
% - GNB single-voxel classifier
%   [accuracyMap,pvalueMap] = computeInformationMap(examples,labels,labelsGroup,'gnb');
%
% - GNB searchlight classifier (searchmight high-speed version)
%   [accuracyMap,pvalueMap] = computeInformationMap(examples,labels,labelsGroup,'gnb_searchmight','searchlight',voxelsToNeighbours,numberOfNeighbours);
%
% - LDA-shrinkage searchlight classifier
%   [accuracyMap,pvalueMap] = computeInformationMap(examples,labels,labelsGroup,'lda_shrinkage','searchlight',voxelsToNeighbours,numberOfNeighbours);
%
% - GNB searchlight classifier (searchmight, permutation test version)
%   [accuracyMap,pvalueMap] = computeInformationMap(examples,labels,labelsGroup,'gnb_searchmight','searchlight',voxelsToNeighbours,numberOfNeighbours,'testToUse','accuracyOneSided_permutation',1000);
%
% Notes:
% - to threshold a p-value map using a clone of Tom Nichols' FDR function FDR.m:
%
% q = 0.01;
% [thresholdID,thresholdN] = computeFDR(pvalueMap,q);  % produces two thresholds, use pN;
% binaryMap = pvalueMap <= thresholdN; % thresholded map
%
% History
% 2010 Apr 29 - fpereira@princeton.edu - thanks to Leila Reddy's courageous beta testing, almost ready
% 2010 Apr 18 - fpereira@princeton.edu - time for a beta!
% 2009 Jul 31 - fpereira@princeton.edu - everything in place, write companion paper
% 2009 Mar 14 - fpereira@princeton.edu - created from previous code
%

function [accuracyMap,pvalueMap,hitMap,falseMap,optionalReturns] = computeInformationMap(varargin)

%% process parameters
this = 'computeInformationMap';
examples    = varargin{1};
labels      = varargin{2};
groupLabels = varargin{3};
classifier  = varargin{4}

testToUse   = 'accuracyOneSided_analytical'; nPermutations = 0;
measure     = 'error';
useSearchlight = 0;
usePriors = 0; % assume that classes are balanced
storeSearchlightCovarianceMatrices = 0;
storeLocalConfusionMatrices = 0;
storeCorrectMatrix = 0;
storeAverageClassMeansStdvs = 0;
computePairmaps = 0;
optionalReturns = [];
addRegressorFeatures = 0;
seed = 1685;
classifierParameters = {};
groupLabelsOriginal = groupLabels;
% Only implemented for searchmight GNB
hitMap = [];
falseMap = [];

if nargin > 4
  % there are additional arguments to process
  % (TODO: test that each <argName> is followed by the right arg)

  idx = 5;
  while idx <= nargin
    argName = varargin{idx}; idx = idx + 1;

    switch argName      
     case {'testToUse'}
      testToUse = varargin{idx}; idx = idx + 1;
      switch testToUse
       case {'accuracyOneSided_permutation'}
        nPermutations = varargin{idx}; idx = idx + 1;
       case {'accuracyOneSided_analytical'}
        % all set
       otherwise
        fprintf('%s: test %s is not supported\n',this,testToUse);return;
      end
     
     case {'errorMeasure'}
      measure = varargin{idx}; idx = idx + 1;
     
     case {'searchlight'}
      % neighbour information given, two arguments intead
      voxelsToNeighbours = varargin{idx}; idx = idx + 1;
      numberOfNeighbours = varargin{idx}; idx = idx + 1;
      useSearchlight = 1;

     case {'usePriors'}
      usePriors = 1;
      
     case {'storeSearchlightCovarianceMatrices'}
      storeSearchlightCovarianceMatrices = 1;
      
     case {'storeLocalConfusionMatrices'}
      storeLocalConfusionMatrices = 1;
      
     case {'storeCorrectMatrix'}
      storeCorrectMatrix = 1;
    
     case {'storeAverageClassMeansStdvs'}
      storeAverageClassMeansStdvs = 1;
      
     case {'computePairmaps'}
      computePairmaps = 1;

     case {'addRegressorFeatures'}
      addRegressorFeatures = 1;
      regressorFeatures = varargin{idx}; idx = idx + 1;

     case {'classifierParameters'}
      % a cell array with all the parameters we would give the classifier code
      classifierParameters = varargin{idx}; idx = idx + 1;
      
     case {'groupLabelsOriginal'}
      % in case we are using even/odd group labels, we need
      % to provide labels to subdivide even/odd folds
      groupLabelsOriginal = varargin{idx}; idx = idx + 1;
      
     case {'seed'}
      seed = varargin{idx}; idx = idx + 1;
    end

  end
end
rand('state',seed); randn('state',seed);

%% argument sanity checking

[n,m] = size(examples);
labelValues = unique(labels); nClasses = length(labelValues);
groupValues = unique(groupLabels); nGroups = length(groupValues);

if nClasses < 2; fprintf('%s: you need to have at least two classes\n',this); return; end
if nGroups < 2;  fprintf('%s: you need to have at least two example groups\n',this); return; end
if nClasses == 2; computePairmaps = 0; end

% convert labels to 1:#classes, as some code depends on that
newLabels = zeros(n,1);
for c = 1:nClasses
  indices = find(labels == labelValues(c));
  newLabels(indices) = c;
end
labels = newLabels; labelValues = 1:nClasses;

for c = 1:nClasses
  maskClass{c}    = (labels == labelValues(c));
  indicesClass{c} = find(maskClass{c});
  nPerClass(c)    = length(indicesClass{c});
end

for g = 1:nGroups
  maskGroup{g}   = (groupLabels == groupValues(g));
  indicesGroup{g} = find(maskGroup{g});
  nPerGroup(g)    = length(indicesGroup{g});
end

switch classifier
 case {'lda','lda_ridge','lda_shrinkage','qda_shrinkage','svm_linear','svm_quadratic','svm_rbf','svm_sigmoid','lr_L1','lr_L2','gnb_searchmight','knn'}
  if ~useSearchlight
    fprintf('%s: in order to use %s you need to specify searchlight\n',this,classifier); return
  end
 case {'gaussian','gaussianPooled','nbayes','nbayesPooled','gnb','gnb_pooled','gnbtest'}
  % all peachy
 otherwise
  fprintf('%s: classifier +%s+ is not supported, argh\n',this,classifier);return;
end

%% some classifiers need the classifier name in the parameters

lcp = length(classifierParameters);
switch classifier
 case {'svm_linear','svm_quadratic','svm_rbf','svm_sigmoid'}
  kernel = classifier(5:end);
  classifierParameters = {'kernel',kernel,classifierParameters{:}};
  lcp = lcp + 2;
 otherwise
end

if useSearchlight
  % check that the neighbour information jigs with the number of voxels
  nNeighbours = size(voxelsToNeighbours,2);
  if (size(voxelsToNeighbours,1) ~= m) || (size(numberOfNeighbours,1) ~= m)
    fprintf('%s: voxelsToNeighbours/numberOfNeighbours do not match #voxels\n',this); return;
  end
end

%
% call the mapping code by building a string containing the command line and then evaluating it
%

% used in both normal and pairmaps

if ~computePairmaps
  % this is the normal case, maps for all pairs of classes are at the bottom 
  
  switch classifier
   case {'gnb_searchmight'}
    % the accuracy map is computed later, as the same function produces the map and permutations
    % (to keep things neat, if analytical p-values are requested it's also done there)
   otherwise
    % complication so that we can have a single command line regardless of the kind of parameters used
    extraParameters = '';
    if useSearchlight
      extraParameters = sprintf('%s,''searchlight'',voxelsToNeighbours,numberOfNeighbours',extraParameters);
    end
    if storeSearchlightCovarianceMatrices
      extraParameters = sprintf('%s,''storeSearchlightCovarianceMatrices''',extraParameters);
    end
    if storeLocalConfusionMatrices
      extraParameters = sprintf('%s,''storeLocalConfusionMatrices''',extraParameters);
    end
    if storeCorrectMatrix
      extraParameters = sprintf('%s,''storeCorrectMatrix''',extraParameters);
    end
    if storeAverageClassMeansStdvs
      extraParameters = sprintf('%s,''storeAverageClassMeansStdvs''',extraParameters);
    end
    if addRegressorFeatures
      %regressorFeatures = varargin{idx}; idx = idx + 1;
      extraParameters = sprintf('%s,''addRegressorFeatures'',regressorFeatures',extraParameters);
    end
    
    if usePriors
      extraParameters = sprintf('%s,''usePriors''',extraParameters);
    end
    
    cmd = sprintf('[sortedFeatures,sortedError,discard,optionalReturns] = computeLocalClassifiers(examples,labels,classifier,''groupLabels'',groupLabels,''classifierParameters'',classifierParameters,''groupLabelsOriginal'',groupLabelsOriginal,''errorMeasure'',''averageRank''');
    if isequal(extraParameters,''); cmd = [cmd,');']; else cmd = sprintf('%s%s);',cmd,extraParameters); end
    eval(cmd);
    
    if 0
      if useSearchlight
        if storeSearchlightCovarianceMatrices
          [sortedFeatures,sortedError] = summarizeELE_rankByNestedCV_g(examples,labels,classifier,'groupLabels',groupLabels,'searchlight',voxelsToNeighbours,numberOfNeighbours,'storeSearchlightCovarianceMatrices');
        else
          [sortedFeatures,sortedError] = summarizeELE_rankByNestedCV_g(examples,labels,classifier,'groupLabels',groupLabels,'searchlight',voxelsToNeighbours,numberOfNeighbours);
        end
      else
        [sortedFeatures,sortedError] = summarizeELE_rankByNestedCV_g(examples,labels,classifier,'groupLabels',groupLabels);
      end
    end
    
    % reorder into original feature order
    [discard,neworder] = sort(sortedFeatures);
    accuracyMap = 1 - sortedError(neworder); % code returns error
  
  end
  
  %
  % compute p-values
  %

  chance = 1/nClasses;
  
  switch testToUse
   case {'accuracyOneSided_analytical'}
    switch classifier 
     case {'gnb_searchmight'}
      % to keep things parallel with the permutation test code, for this classifier the maps are obtained here
      [discard,newOrder] = sort(groupLabels); % sort so that the groups come in order (makes code simpler)
      neighbourMax    = size(voxelsToNeighbours,2);
      neighbourRadius = round((neighbourMax+1)^(1/3)-1)/2;
      
      [accuracyMap,discard,pvalueMap,hitMap,falseMap] = searchmightGNB(examples(newOrder,:)',labels(newOrder)',groupLabels(newOrder)',neighbourRadius,voxelsToNeighbours',numberOfNeighbours',0); clear discard;
      %accuracyMap = accuracyMap'; pvalueMap = pvalueMap';
     otherwise
      % we already have the accuracyMap
    end
      % convert accuracy to count of #examples out of n that were labelled correctly
    countMap  = round(accuracyMap*n);  
    
    % P(X>=observed|H0 is true) = 1 - P(X<observed|H0 is true)
    pvalueMap = 1-(binocdf(countMap,n,chance)-binopdf(countMap,n,chance));
    
   case {'accuracyOneSided_permutation'}
 
    switch classifier
     case {'gnb_searchmight'}
      % requires a completely different command line
      [discard,newOrder] = sort(groupLabels); % sort so that the groups come in order (makes code simpler)
      neighbourMax    = size(voxelsToNeighbours,2);
      neighbourRadius = round((neighbourMax+1)^(1/3)-1)/2;
      
      [accuracyMap,discard,pvalueMap,hitMap,hitMap] = searchmightGNB(examples(newOrder,:)',labels(newOrder)',groupLabels(newOrder)',neighbourRadius,voxelsToNeighbours',numberOfNeighbours',nPermutations); clear discard;
      accuracyMap = accuracyMap'; pvalueMap = pvalueMap';
     
     otherwise
      labelsPermuted  = zeros(size(labels));
      aboveOrEqualMap = zeros(1,m); % count how often accuracy on permuted data goes over observed
      fprintf('computing permutation test p-values\n');
      
      % prepare commandline
      cmd = sprintf('[sortedFeatures,sortedError,discard,optionalReturns] = computeLocalClassifiers(examples,labelsPermuted,classifier,''groupLabels'',groupLabels','useSilence');
      if isequal(extraParameters,''); cmd = [cmd,');']; else; cmd = sprintf('%s%s);',cmd,extraParameters); end
    
      % compute
      countEqual = 0; % # of permutations equal to original labelling
    
      for p = 1:nPermutations
        if (rem(p,10) == 1) tic; end
        
        % permute labels within each group
        for g = 1:nGroups
          labelsPermuted(indicesGroup{g}) = labels(indicesGroup{g}(randperm(nPerGroup(g))));
        end

        if labelsPermuted == labels; countEqual = countEqual + 1; end
      
        % call the prepared commandline
        eval(cmd);
      
        %    if useSearchlight
        %      [sortedFeatures,sortedError] = computeLocalClassifiers(examples,labelsPermuted,classifier,'groupLabels',groupLabels,'searchlight',voxelsToNeighbours,numberOfNeighbours,'useSilence');
        %    else
        %      [sortedFeatures,sortedError] = computeLocalClassifiers(examples,labelsPermuted,classifier,'groupLabels',groupLabels,'useSilence');
        %    end
        
        % reorder into original feature order
        [discard,neworder] = sort(sortedFeatures);
        accuracyMapPermuted = 1 - sortedError(neworder); % code returns error
      
        if (rem(p,10) == 1) t=toc; fprintf('\tcomputed permutation %d/%d: estimated time to completion: %s seconds\n',p,nPermutations,num2str(t*(nPermutations-p)));end
      
        % tally whether accuracy on permuted data went over observed
        indicesHigher = find(accuracyMapPermuted >= accuracyMap);
        aboveOrEqualMap(indicesHigher) = aboveOrEqualMap(indicesHigher) + 1;
      end; % end of for over permutations

      % compute p-values from tally
      % the value is the fraction of the number of permutations where the accuracy was >= observed accuracy
      pvalueMap = aboveOrEqualMap / nPermutations;
    
    end; % of switch on classifier 
  end; % of switch on test to use


else
  %
  % computation of maps for each pair of classes
  % 

  nPairs = nClasses*(nClasses-1)/2;
  
  [n,m] = size(examples);
  labelValues = unique(labels); nClasses = length(labelValues);
  groupValues = unique(groupLabels); nGroups = length(groupValues);
  
  accuracyMap     = zeros(nPairs,m);
  pvalueMap       = zeros(nPairs,m);
  optionalReturns = cell(nPairs,1);
  
  p = 1;
  for c1 = 1:(nClasses-1)
    indicesC1 = indicesClass{c1};
    for c2 = (c1+1):nClasses
      indicesC2 = indicesClass{c2};
      
      % call this for the examples in this pair of classes
      
      indicesPair = [indicesC1;indicesC2];
      
      cmd = sprintf('[accuracyMap(p,:),pvalueMap(p,:),optionalReturns{p}]=computeInformationMap(examples(indicesPair,:),labels(indicesPair),groupLabels(indicesPair),classifier,''classifierParameters'',classifierParameters');
      for i = 5:nargin
        cmd = sprintf('%s,varargin{%d}',cmd,i);
      end
      cmd = sprintf('%s);',cmd);
      eval(cmd);
      
      p = p + 1;
    end
  end  
end


function [] = testThis() 

%% simple two class dataset with a 3D mask

dimx = 11; dimy = 11; dimz = 12;

mask = zeros(dimx,dimy,dimz);
xrange = 3:7; yrange = 3:7;
for iz = 2:(dimz-1); mask(xrange,yrange,iz) = 1; end
meta = createMetaFromMask(mask);
nVoxels = length(meta.indicesIn3D);

% define a few voxels with different means in different conditions
if 0
  nClasses = 2;

  voxelMeans3D = zeros(dimx,dimy,dimz,nClasses);
  voxelMeans3D(5,5,2,2) = 2;
  voxelMeans3D(5,5,2,2) = 2;
  voxelMeans3D(5,5,3,2) = 2;
  voxelMeans3D(5,5,3,2) = 2;
  voxelMeans3D(5,5,6,2) = 1;
  voxelMeans3D(5,5,6,2) = 1;
  voxelMeans3D(5,5,7,2) = 1;
  voxelMeans3D(5,5,7,2) = 1;
else
  nClasses = 3;
  voxelMeans3D = zeros(dimx,dimy,dimz,nClasses);

  voxelMeans3D(5,5,2,2) = 2;
  voxelMeans3D(5,5,2,2) = 2;
  voxelMeans3D(5,5,3,2) = 2;
  voxelMeans3D(5,5,3,2) = 2;
  voxelMeans3D(5,5,6,2) = 1;
  voxelMeans3D(5,5,6,2) = 1;
  voxelMeans3D(5,5,7,2) = 1;
  voxelMeans3D(5,5,7,2) = 1;
  
  voxelMeans3D(5,5,2,3) = 4;
  voxelMeans3D(5,5,2,3) = 4;
  voxelMeans3D(5,5,3,3) = 4;
  voxelMeans3D(5,5,3,3) = 4;
  voxelMeans3D(5,5,6,3) = 2;
  voxelMeans3D(5,5,6,3) = 2;
  voxelMeans3D(5,5,7,3) = 2;
  voxelMeans3D(5,5,7,3) = 2;
end
  
for ic = 1:nClasses
  tmp = voxelMeans3D(:,:,:,ic);
  voxelMeans{ic} = tmp(meta.indicesIn3D)';
end
  
nGroups = 4;
nPerClassPerGroup = 25; npcpg = nPerClassPerGroup;

n = nGroups * nPerClassPerGroup * 2;
m = nVoxels;

examples    = zeros(n,m);
labels      = zeros(n,1);
labelsGroup = zeros(n,1);

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

imagesc([examples,labels])

classifier = 'gnb_pooled';
[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier);

classifier = 'gnb_pooled';
classifier = 'lda_shrinkage';
classifier = 'svm_linear';
classifier = 'svm_rbf';
classifier = 'svm_quadratic';
classifier = 'svm_sigmoid';
[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours);


classifier = 'svm_linear';
classifierParameters = {'lambda',1,'gamma',1};
%classifierParameters = {'lambda','crossValidation','gamma',1};
classifierParameters = {'lambda','crossValidation'};
classifier = 'svm_rbf';
classifier = 'svm_quadratic';
classifier = 'svm_sigmoid';
classifierParameters = {'lambda','crossValidation','gamma','crossValidation'};

[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours,'classifierParameters',classifierParameters);

labelsGroupEvenOdd = rem(labelsGroup,2) + 1;
[am,pm] = computeInformationMap(examples,labels,labelsGroupEvenOdd,classifier,'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours,'classifierParameters',classifierParameters,'groupLabelsOriginal',labelsGroup);

[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours,'classifierParameters',classifierParameters,'computePairmaps');

clf;
slices = 1:8; nSlices = length(slices);
volume = repmat(NaN,[dimx dimy dimz]);
volume(meta.indicesIn3D) = am;

nrows = 3; ncols = 3;
for iz = 1:nSlices
  z = slices(iz);
  subplot(nrows,ncols,iz);
  imagesc(volume(:,:,z),[0 1]); axis square;
end

%%
%% KNN requires examples created by corrupting a per class pattern
%%

nClasses = 3;
templateShared = randn(4,4);
for ic = 1:nClasses; templateClass{ic} = randn(4,4); end
dimz = nClasses*2 + 1; % class specific patterns end up in slices separated by 1
dimx = 8; dimy = 8; % needs enough room to paste
xrange = 1:4; yrange = 1:4;

mask = ones(dimx,dimy,dimz);
meta = createMetaFromMask(mask);
nVoxels = length(meta.indicesIn3D);

voxelMeans3D = zeros(dimx,dimy,dimz,nClasses);
iz = 3;
for ic = 1:nClasses
  voxelMeans3D(xrange,yrange,1,ic)  = templateShared;
  voxelMeans3D(xrange,yrange,iz,ic) = templateClass{ic};
  iz = iz + 2;
end
  
for ic = 1:nClasses
  tmp = voxelMeans3D(:,:,:,ic);
  voxelMeans{ic} = tmp(meta.indicesIn3D)';
end
  
nGroups = 4;
nPerClassPerGroup = 25; npcpg = nPerClassPerGroup;

n = nGroups * nPerClassPerGroup * 2;
m = nVoxels;

if 0
  % plot for debugging
  clf; nrows = nClasses; ncols = dimz;
  idx = 1;
  for ic = 1:nClasses
    for iz = 1:dimz
      subplot(nrows,ncols,idx);
      imagesc(voxelMeans3D(:,:,iz,ic)); axis square;
      idx = idx + 1;
    end
  end
  pause
end

examples    = zeros(n,m);
labels      = zeros(n,1);
labelsGroup = zeros(n,1);

idx = 1;
for ig = 1:nGroups
  for ic = 1:nClasses
    range = idx:(idx+npcpg-1);
    
    examples(range,:) = randn(npcpg,m) * 0.25;
    examples(range,:) = examples(range,:) + repmat(voxelMeans{ic},npcpg,1);
    
    labels(range)      = ic;
    labelsGroup(range) = ig;
    idx = idx + npcpg;
  end
end

labels = labels + 1;


classifier = 'knn'; classifierParameters = {'k',0,'similarityMeasure','correlation'};
classifier = 'knn'; classifierParameters = {'k',1,'similarityMeasure','correlation'};

[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours,'classifierParameters',classifierParameters);



classifier = 'svm_linear'; classifierParameters = {'lambda',1,'gamma',1};
[am,pm] = computeInformationMap(examples,labels,labelsGroup,classifier,'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours,'classifierParameters',classifierParameters);
