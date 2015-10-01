%
% Draft of a function to compute many maps from a single dataset, using nice names
%
% In:
% - examples
% - labels
% - labelsGroup
% - optional:
%
% Out:
% - array of accuracy maps
% - array of p-value maps
%
% History
% 2009 April 3 - fpereira@princeton.edu - created
%

function [accuracyMaps,pvalueMaps,voxelRankings] = computeManyMapsFromOneDataset(varargin)

%% process parameters

this = 'computeManyMapsFromOneDataset';
examples    = varargin{1}; [nExamples,nVoxels] = size(examples);
labels      = varargin{2};
labelsGroup = varargin{3};
meta        = varargin{4};

% derive a few things that are useful in setting the defaults
tmp = unique(labels); nClasses = length(tmp);
labelsNew = zeros(size(labels)); for c = 1:nClasses; indices = find(labels == tmp(c)); labelsNew(indices) = c; end;
labels = labelsNew; labelValues = 1:nClasses;

groupValues = unique(labelsGroup); nGroups = length(groupValues);

if nClasses < 2; fprintf('%s: you need to have at least two classes\n',this); return; end
if nGroups < 2;  fprintf('%s: you need to have at least two example groups\n',this); return; end

% set defaults
testType = 'accuracyOneSided_analytical';
storeLocalInformation = 0;
computeRankings = 0; voxelRankings = [];

if nClasses == 2
  maps = {'voxelwiseGNB','searchlightGNB','searchlightLDA_shrinkage','searchlightLDA_ridge','svm_linear','svm_quadratic','svm_rbf'};
else
  % skip the SVM maps
  maps = {'voxelwiseGNB','searchlightGNB','searchlightLDA_shrinkage','searchlightLDA_ridge'};
end

% now get optional arguments

if nargin > 4
  % there are additional arguments to process

  idx = 5;
  while idx <= nargin
    argName = varargin{idx}; idx = idx + 1;

    switch argName  
     case {'maps'}
      maps = varargin{idx}; idx = idx + 1;
    end

  end
end

% and check
nMaps = length(maps);

for m = 1:nMaps
  
  switch maps{m}
   case {'voxelwiseGNB','voxelwiseGNBsmoothed','searchlightGNB','searchlightLDA_shrinkage','searchlightLDA_ridge','searchlightSVM_linear','searchlightSVM_quadratic','searchlightSVM_rbf'}
    % OK
   otherwise
    fprintf('%s: map type %s is not supported\n',this,maps{m}); return;
  end
end

%
% Compute
%

accuracyMaps  = zeros(nMaps,nVoxels);
pvalueMaps    = zeros(nMaps,nVoxels);
voxelRankings = []; if computeRankings; voxelRankings = zeros(nMaps,nVoxels); end
  
for m = 1:nMaps
  
  map = maps{m}; fprintf('%s: computing %s map\n',this,map);
  optionalReturns{m} = [];
  
  % convert nice map name to classifier name in the function
  switch map
   case {'voxelwiseGNB','voxelwiseGNBsmoothed','searchlightGNB'}
    classifier = 'nbayesPooled';
   case {'searchlightGNB','searchlightGNBsmoothed'}
    classifier = 'nbayesPooled';
   case {'searchlightLDA'}
    classifier = 'lda';
   case {'searchlightLDA_ridge'}
    classifier = 'lda_ridge';
   case {'searchlightLDA_shrinkage'}
    classifier = 'lda_shrinkage';
   case {'searchlightSVM_linear'}
    classifier = 'svm_linear';
   case {'searchlightSVM_quadratic'}
    classifier = 'svm_quadratic';
   case {'searchlightSVM_rbf'}
    classifier = 'svm_rbf';
   otherwise
    fprintf('error: should not be here...\n');return
  end
  
  switch map
    % single voxel classifiers
   case {'voxelwiseGNB','voxelwiseGNBsmoothed'}

    switch map
     case {'voxelwiseGNBsmoothed'}
      [examplesHere] = simplesmooth3D(examples,meta.voxelsToNeighbours,meta.numberOfNeighbours,'smoothingCriterion','box',{0});
     otherwise
      examplesHere = examples;
    end
    
    switch testType
     case {'accuracyOneSided_analytical'}
      [accuracyMaps(m,:),pvalueMaps(m,:),optionalReturns{m}] = computeInformationMap(examplesHere,labels,labelsGroup,classifier);
     case {'accuracyOneSided_permutation'}
      [accuracyMaps(m,:),pvalueMaps(m,:)] = computeInformationMap(examplesHere,labels,labelsGroup,'nbayesPooled','testToUse','accuracyOneSided_permutation',nPermutations);
    end
    clear examplesHere;

    % searchlight classifiers
    
   case {'searchlightGNB','searchlightGNBsmoothed','searchlightLDA','searchlightLDA_ridge','searchlightLDA_shrinkage','searchlightSVM_linear','searchlightSVM_quadratic','searchlightSVM_rbf'}
    
    switch map
     case {'searchlightGNBsmoothed'}
      % smooth data classifiers
      [examplesHere] = simplesmooth3D(examples,meta.voxelsToNeighbours,meta.numberOfNeighbours,'smoothingCriterion','box',{0});
     otherwise
      examplesHere = examples;
    end

    switch testType
     case {'accuracyOneSided_analytical'}
      if storeLocalInformation
        if nClasses == 2
          [accuracyMaps(m,:),pvalueMaps(m,:),optionalReturns{m}] = computeInformationMap(examplesHere,labels,labelsGroup,classifier,...
                                                            'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours,'storeSearchlightCovarianceMatrices');
        else
          [accuracyMaps(m,:),pvalueMaps(m,:),optionalReturns{m}] = computeInformationMap(examplesHere,labels,labelsGroup,classifier,...
                                                            'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours,'storeLocalConfusionMatrices');         
        end
      else
        [accuracyMaps(m,:),pvalueMaps(m,:),optionalReturns{m}] = computeInformationMap(examplesHere,labels,labelsGroup,classifier,'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours);
      end
        
     case {'accuracyOneSided_permutation'}                            
      [accuracyMaps(m,:),pvalueMaps(m,:)] = computeInformationMap(examplesHere,labels,labelsGroup,'nbayesPooled',...
                                             'searchlight',meta.voxelsToNeighbours,meta.numberOfNeighbours,...
                                             'testToUse','accuracyOneSided_permutation',nPermutations);
    end; % of swith on test
    clear examplesHere;
    
  end; % of switch on map

  if computeRankings
    % compute ranking (just sort the accuracy)
    [discard,voxelRankings(m,:)] = sort(accuracyMaps(m,:),'descend');
  end
    
end; % of for loop over maps
