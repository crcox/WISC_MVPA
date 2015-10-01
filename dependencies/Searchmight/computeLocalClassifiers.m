% Ranks features by using each one as a gaussian bayesian
% classifier in cross-validation within the training set.
%
% The code defaults to leaving one example of each class out in every fold,
% but it can also operate by indicating a group membership for each example,
% and leaving out each group of examples in turn.
%%
% In:
% - examples (#examples x #features array)
% - labels   (#examples x 1 vector)
% - classifier to use:
%   - either 'gaussian' or 'gaussianPooled' (assumes same variance for both classes)
% - optional arguments, in any order
%
%   - 'classifierParameters',<classifier parameters cell array>
%      if the classifier requires parameters, computeInformationMap will handle parameter
%      processing and checking and pass them in this cell array
%
%   - 'errorMeasure',<'error'|'averageRank'> - error measure to use 
%     (defaults to 'error', 'averageRank' is the normalized average
%      rank in a multiclass situation:
%      for every example, the classifier outputs a ranking of possible classes
%      the measure is the position in the rank of the correct
%      label, normalized rank is just this measure normalized to [0,1]
%      This is equivalent to error, in a 2-class situation.
%     )
%
%   - 'groupLabels', <#examples x 1>  - vector with each example's group membership for crossvalidation
%     (defaults to leaving one example of each class out)
%     (e.g. for each example, have an indicator of what run it comes from and the
%      the code will compute leave-one-run-out-crossvalidation)
%
%   - 'searchlight',<voxelsToNeighbours>,<numberOfNeighbours> - 
%      compute the error measure for a searchlight centred on the voxels, with searchlight
%      details being given by the neighbourhood.
%
%   - 'errorMeasure',<'error'|'averageRank'> - error measure to use 
%     (defaults to 'error', 'averageRank' is the normalized average
%      rank in a multiclass situation:
%      for every example, the classifier outputs a ranking of possible classes
%      the measure is the position in the rank of the correct
%      label, normalized rank is just this measure normalized to [0,1]
%      This is equivalent to error, in a 2-class situation.
%     )
%
%     DO NOT USE THESE UNLESS YOU KNOW WHAT YOU ARE DOING
%
%   - 'usePriors' - by default, the code assumes the classes are balanced and thus their prior
%                   probabilities are the same and need not be considered. This switch brings
%                   the priors into the calculation.
%
%   - 'useWholeTrainingSet' - if specified, do not cross-validate and
%                          use the entire training set for both train and test

%
% Out:
% - sortedFeatures  - features sorted from best to worst
% - sortedMeasure   - feature errorMeasure sorted from best to worst (matches sortedFeatures)
% - measurePerClass - the measure for each example, separated by classes
%
% Notes:
% - assumes a balanced dataset, i.e. an equal number of examples of each class
% - extensive test code at the end of the function
%
% Examples:
% - [sortedFeatures,sortedMeasure] = summarizeELE_rankByNestedCV_g( examples,labels,'gaussianPooled')
% - [sortedFeatures,sortedMeasure] = summarizeELE_rankByNestedCV_g( examples,labels,'gaussianPooled','groupLabels',glabels)
% - [sortedFeatures,sortedMeasure] = summarizeELE_rankByNestedCV_g( examples,labels,'gaussianPooled','groupLabels',glabels,'searchlight',voxelsToNeighbours,numberOfNeighbours)
%
% History
% - 2009 April 1 - fpereira@princeton.edu - created from previous code, together with computeInformationMap
%

function [sortedFeatures,sortedMeasure,measurePerClass,optionalReturns] = computeLocalClassifiers( varargin )

this = 'computeLocalClassifiers';

%% process arguments and figure out/test a few things

if nargin < 3; eval(sprintf('help %s',this));return; end

examples   = varargin{1};
labels     = varargin{2};
classifier = varargin{3};
classifierParameters = {};

%maps = {'voxelwiseGNB','voxelwiseGNBsmoothed','searchlightGNB','searchlightLDA_shrinkage','searchlightLDA_ridge'};
%maps = {'voxelwiseGNB','searchlightGNB','searchlightLDA_shrinkage'};
%maps = {'voxelwiseGNB','searchlightGNB','searchlightLDA_shrinkage','searchlightSVM_linear','searchlightSVM_quadratic','searchlightSVM_rbf'};
%maps = {'searchlightLR'};
%nmaps = length(maps);

switch classifier
 case {'gaussian','gaussianPooled','nbayes','nbayesPooled','gnb','gnb_pooled','gnbtest','lda','lda_ridge','lda_shrinkage','qda_shrinkage','svm_linear','svm_quadratic','svm_sigmoid','svm_rbf','lr_L1','lr_L2','searchmightGNB','knn'}
  % OK
 otherwise
  fprintf('%s: classifier %s is not supported\n',this,classifier);pause;return;
end

sortedLabelValues = unique(labels); nLabels = length(sortedLabelValues);
[nExamples,nFeatures] = size(examples);

% turn labels into the 1:#classes range
labelsNew = zeros(size(labels));
for c = 1:nLabels; indices = find(labels == sortedLabelValues(c)); labelsNew(indices) = c; end; clear indices;
labels = labelsNew; sortedLabelValues = 1:nLabels;

measure     = 'error';
groupLabels = []; groupLabelsOriginal = [];
useSearchlight = 0;
useWholeTrainingSet = 0;
usePriors = 0; % assume that classes are balanced
addRegressorFeatures = 0;
useCount  = 1;
useSilence = 0;
optionalReturns.searchlightCovarianceMatrices = []; storeSearchlightCovarianceMatrices = 0;
optionalReturns.localConfusionMatrices = [];        storeLocalConfusionMatrices = 0;
optionalReturns.correctMatrix = [];                 storeCorrectMatrix = 0;
optionalReturns.averageClassMeans = []; optionalReturns.averageClassStdvs = []; storeAverageClassMeansStdvs = 0;
classifierParameters = {};

if nargin > 3
  % there are additional arguments to process
  % (TODO: test that each <argName> is followed by the right arg)

  idx = 4;
  while idx <= nargin
    argName = varargin{idx}; idx = idx + 1;

    switch argName      
     case {'errorMeasure'}
      measure = varargin{idx}; idx = idx + 1;
     
     case {'groupLabels'}
      groupLabels = varargin{idx}; idx = idx + 1;
     
     case {'searchlight'}
      if isinteger(varargin{idx})
	% we're given radius, compute neighbours explicitly
	radius     = varargin{idx}; idx = idx + 1;
	colToCoord = varargin{idx}; idx = idx + 1;
	computeNeighbourInformation = 1;
      else
	% neighbour information given, two arguments intead
	voxelsToNeighbours = varargin{idx}; idx = idx + 1;
	numberOfNeighbours = varargin{idx}; idx = idx + 1;
	computeNeighbourInformation = 0;
      end
      useSearchlight = 1;

     case {'usePriors'}
      usePriors = 1;
      
     case {'useWholeTrainingSet'}
      useWholeTrainingSet = 1; 
      
     case {'storeSearchlightCovarianceMatrices'}
      % stores the average local covariance matrix (average over folds)
      storeSearchlightCovarianceMatrices = 1;
      optionalReturns.searchlightCovarianceMatrices = zeros(27,27,nFeatures); % some will be smaller than 27x27
      optionalReturns.searchlightMeans = zeros(nFeatures,27,nLabels);
      if ~useSearchlight
        fprintf('%s: storeSearchlightCovarianceMatrices can only be used with LDA classifiers\n',this);return;
      end

     case {'storeLocalConfusionMatrices'}
      % stores the local confusion matrix (searchlight or not)
      storeLocalConfusionMatrices = 1;
      optionalReturns.localConfusionMatrices = zeros(nLabels,nLabels,nFeatures);
     
     case {'storeCorrectMatrix'}
      % this is the #examples x #voxels matrix of binary values
      storeCorrectMatrix = 1;
      optionalReturns.correctmatrix = zeros(nExamples,nFeatures);
   
     case {'storeAverageClassMeansStdvs'}
      storeAverageClassMeansStdvs = 1;
      optionalReturns.averageClassMeans = zeros(nLabels,nFeatures);
      optionalReturns.averageClassStdvs = zeros(nLabels,nFeatures);
      
     case {'addRegressorFeatures'}
      addRegressorFeatures = 1;
      regressorFeatures = varargin{idx}; idx = idx + 1;

     case {'classifierParameters'}
      % a cell array with all the parameters we would give the classifier code
      classifierParameters = varargin{idx}; idx = idx + 1;
      
     case {'useSilence'}
      % cuts down most messages
      useSilence = 1;
    
     case {'groupLabelsOriginal'}
      groupLabelsOriginal = varargin{idx}; idx = idx + 1;    
    end

  end
end

if isempty(groupLabelsOriginal) & ~isempty(groupLabels)
  groupLabelsOriginal = groupLabels;
end

if addRegressorFeatures
  % add mock neighbours to every voxel corresponding to all the new regressors
  [tmp,nRegressorFeatures] = size(regressorFeatures);
  if tmp ~= nExamples; fprintf('%s: regressor features must have same #rows as data\n',this); return; end
    
  %numberOfNeighbours = numberOfNeighbours + 1;
  maxNeighbours = size(voxelsToNeighbours,2);
  
  voxelsToNeighbours = [voxelsToNeighbours,zeros(nFeatures,nRegressorFeatures)];
    
  newNeighbours = nFeatures + (1:nRegressorFeatures);
  for v = 1:nFeatures
    nrange = numberOfNeighbours(v) + (1:nRegressorFeatures);
    voxelsToNeighbours(v,nrange) = newNeighbours;
  end
  numberOfNeighbours = numberOfNeighbours + nRegressorFeatures;;
  
  examples = [examples,regressorFeatures];
end


%% extra searchlight processing

if useSearchlight
  
  if isempty(groupLabels)
    fprintf('%s: in order to use searchlight you must specify group labels\n',this); return
  end

  if computeNeighbourInformation
    fprintf('%s: computing neighbour information, consider precomputing this\n',this); return      
    [voxelsToNeighbours,numberOfNeighbours] = neighboursWithinRadius(colToCoord,radius);
  end
else
  if isequal(classifier,'lda') | isequal(classifier,'lda_ridge') | isequal(classifier,'lda_shrinkage')| isequal(classifier,'qda_shrinkage')
    fprintf('%s: in order to use lda you need to specify searchlight\n',this); return
  end
end

%% classifier specific preprocessing

switch classifier
 case {'svm_linear','svm_quadratic','svm_sigmoid','svm_rbf'}
  %% if using libsvm, the code will order labels by the first time they appear, so
  %% things are simple if, within each group, examples are sorted by class
  %% (ideally this should be fixed directly in the label prediction wrapper code)
  
  if 1
    groups = unique(groupLabels); nGroups = length(groups);
    examplesNew    = zeros(size(examples));
    labelsNew      = zeros(size(labels));
    groupLabelsNew = zeros(size(groupLabels));
  
    for g = 1:nGroups
      indices = find(groupLabels == groups(g));
      [discard,neworder] = sort(labels(indices));
      examplesNew(indices,:)  = examples(indices(neworder),:);
      labelsNew(indices)      = labels(indices(neworder));
      groupLabelsNew(indices) = groupLabels(indices(neworder));
    end
    examples    = examplesNew;    clear examplesNew;
    labels      = labelsNew;      clear labelsNew;
    groupLabels = groupLabelsNew; clear groupLabelsNew;
  end
  
 otherwise
end
  
%% Find the indices of the examples in each class

indicesLabel  = cell(nLabels,1);
nPerLabel     = zeros(nLabels,1);

for l = 1:nLabels  
  % find examples with this label
  label = sortedLabelValues(l); 
  maskLabel{l}    = (labels == label);
  indicesLabel{l} = find(maskLabel{l});
  nPerLabel(l)    = length(indicesLabel{l});
end

% ensure we have a balanced dataset
if sum(diff(nPerLabel))
  unbalanced = 1;
else
  unbalanced = 0;
end

%
% The cross-validation code has two different branches, depending on whether
% group labels were given (leave-1-group-out) or not (leave-1-example-of-each-class-out).
% The latter can be done efficiently with a few tricks, but the code is more complicated.
% Leaving one group out is simpler and is the branch that comes first.
%

% hack to force us to use the group code
if isempty(groupLabels)
  if unbalanced
    % the rest of the code can handle unbalanced if this part is fine
    fprintf('error: there must be the same number of examples in each class\n');pause;return;
  else
    groupLabels = zeros(nExamples,1);
    for l = 1:nLabels
      groupLabels(indicesLabel{l}) = (1:nPerLabel(1))';
    end
  end
end

% The easiest way to use the whole training set is to make all the
% examples be in the same group. Inside the group loop a hack will
% make train and test be the same group.

if useWholeTrainingSet; groupLabels = ones(nExamples,1); end

%
% leave-1-group-out crossvalidation
%
groups = unique(groupLabels); nGroups = length(groups);

% used to store the scores over labels attributed by each feature
% to each example, but it uses too much memory if there are many labels
%  scores = zeros(nExamples,nFeatures,nLabels);

measurePerExample = zeros(nExamples,nFeatures);

if ~useSilence
  if ~useSearchlight
    fprintf('%s: compute leave-1-group-out crossvalidation error measure, voxelwise\n',this);
  else
    fprintf('%s: compute leave-1-group-out crossvalidation error measure, using searchlights\n',this);
  end
end

% some classifiers need initialization of output structures
switch classifier
 case {'svm_linear','svm_quadratic','svm_sigmoid','svm_rbf'}
  optionalReturns.lambdas = zeros(nGroups,nFeatures);
  optionalReturns.gammas  = zeros(nGroups,nFeatures);
 otherwise
end
  
%
% main cross-validation loop
%

for g = 1:nGroups
  group = groups(g); % group to leave out
  if ~useSilence; fprintf('\tprocessing group %d\n',g); end
  
  mask = (groupLabels == group);
  indicesTest  = find( mask); nTest  = length(indicesTest);
  indicesTrain = find(~mask); nTrain = length(indicesTrain);
  
  exampleMeanTrainClass = zeros(nLabels,nFeatures);
  exampleStdvTrainClass = zeros(nLabels,nFeatures);
  
  for l = 1:nLabels
    indicesTestClass{l}  = find(  mask & maskLabel{l} );
    indicesTrainClass{l} = find( ~mask & maskLabel{l} );
  
    exampleMeanTrainClass(l,:) = mean(examples(indicesTrainClass{l},:),1);
    exampleStdvTrainClass(l,:) = std( examples(indicesTrainClass{l},:),0,1);
  end
  
  % make train and test data be the same
  if useWholeTrainingSet
    indicesTrain = indicesTest; nTrain = nTest;
  end

  labelsTest    = labels(indicesTest);    
  labelsTrain   = labels(indicesTrain);

  scoresTest = zeros(nTest,nFeatures,nLabels);

  switch classifier

%%%% GNB classifiers 
    
   case {'gaussian','nbayes','gaussianPooled','nbayesPooled','gnb','gnbtest','gnb_pooled'}
    
    %% compute the class mean and standard deviation for each feature 
    %% inside the training set
    
    if useSearchlight;
      groupMean = cell(nLabels,nFeatures); 
      groupStdv = cell(nLabels,nFeatures); 
    end

    nTrainPerLabel = zeros(1,nLabels);
    for l = 1:nLabels
      label = sortedLabelValues(l);
      
      indices{l}   = find(labelsTrain == label); % training set examples with this label
      groupMean{l} = mean(examples(indicesTrain(indices{l}),:),1);
      nTrainPerLabel(l) = length(indices{l});
      
      switch classifier
       case {'gaussian','nbayes','gnb'}
        % each class has its own standard deviation estimate
        groupStdv{l} = std(examples(indicesTrain(indices{l}),:),0,1);
       case {'gaussianPooled','nbayesPooled','gnb_pooled','lda','lda_ridge','lda_shrinkage'}
        % subtract the mean from these examples, so that we can use
        % them and examples with other labels to compute the standard deviation
        examples(indicesTrain(indices{l}),:) = examples(indicesTrain(indices{l}),:)-repmat(groupMean{l},nTrainPerLabel(l),1);
       case {'gnbtest'}
        % no mean subtraction
      end
    end; % for over labels
    classPriorsTrain = nTrainPerLabel / sum(nTrainPerLabel);
    
    switch classifier
     case {'gaussian','nbayes','gnb'}
      % all done
     case {'gaussianPooled','nbayesPooled','gnb_pooled'}
      % compute standard deviation using all examples after mean subtraction
      tmp = std(examples(indicesTrain,:),0,1);
      for l = 1:nLabels; groupStdv{l} = tmp; end
     case {'gnbtest'}
      % NOTE: working with variance instead, as that is all that is needed
      %tmp = sum(examples(indicesTrain,:).^2,1)/(nTrain-1); % variance
      %tmp = std(examples(indicesTrain,:),0,1);

      % compute it like searchmightGNB.c
      tmp = zeros(1,nFeatures);
      for l = 1:nLabels
        tmp = tmp + sum( examples(indicesTrain(indices{l}),:).^2 - ...
                         2*repmat(groupMean{l},nTrainPerLabel(l),1).*examples(indicesTrain(indices{l}),:) +...
                         repmat(groupMean{l},nTrainPerLabel(l),1).^2, 1);
      end
      tmp = tmp / (nTrain-1);
      
      for l = 1:nLabels; groupStdv{l} = tmp; end
     case {'lda','lda_ridge','lda_shrinkage'}
      % any reasonable way of computing voxelwise parts of covariance matrix?
      
      % Although different matrices have different dimensions, it would be more
      % efficient to store in a conceptual 27x27x #voxels matrix. In practice,
      % it will be 
    end

    if storeAverageClassMeansStdvs
      for l = 1:nLabels
        optionalReturns.averageClassMeans(l,:) = optionalReturns.averageClassMeans(l,:) + groupMean{l};
        optionalReturns.averageClassStdvs(l,:) = optionalReturns.averageClassStdvs(l,:) + groupStdv{l};
      end
    end
    
    % now restore the means, if necessary
    
    for l = 1:nLabels
      label = sortedLabelValues(l);
      
      %indices{l} = find(labelsTrain == label); % training set examples with this label
      
      switch classifier
       case {'gaussian','nbayes','gnb'}
        % nothing to do
       case {'gaussianPooled','nbayesPooled','gnb_pooled'}
        % re-add the class means
        examples(indicesTrain(indices{l}),:) = examples(indicesTrain(indices{l}),:)+repmat(groupMean{l},nTrainPerLabel(l),1);
       case {'gnbtest'}
      end
      
    end; % for over labels
    
    
    %% use the learnt means and standard deviations to compute the
    %% class score for each feature on the test examples

    % class priors train is 1 x #classes
    if usePriors;  matrixOfPriors = log(repmat(classPriorsTrain,nTest,1)); end
    
    for lA = 1:nLabels
      
      if ~useSearchlight
        
        switch classifier
         case {'gaussian','gaussianPooled','nbayes','nbayesPooled','gnb','gnb_pooled'}
          % compute
          % log(prob(class<lA>|feature)), for all classes <lA>
          % assume priors are equal for all classes (given that
          % dataset is balanced)
          
          tmp = examples(indicesTest,:);
          tmp = tmp -  repmat(groupMean{lA},nTest,1);
          tmp = tmp ./ repmat(groupStdv{lA},nTest,1);
          tmp = -1 * (tmp .^2) / 2;
          tmp = tmp - log(sqrt(2*pi)*repmat(groupStdv{lA},nTest,1));

          if usePriors; tmp = tmp + matrixOfPriors(:,lA); end
          
          % use to store scores directly, too much memory in a single structure
          % scores(indicesTest,:,lA) = tmp;  

          % tmp is #test examples x #features
          
          scoresTest(:,:,lA) = tmp; clear tmp;
          
         case {'gnbtest'}
          fprintf('ERROR: not implemented yet\n');pause;
          % needs to be implemented to take into account our keeping variances
        end 
        
      else
        
        switch classifier
         case {'gaussian','gaussianPooled','nbayes','nbayesPooled','gnb','gnb_pooled'}
          
          if ~useSilence; fprintf('\tgroup %d: label %d: searchlight applying feature by feature\n',g,lA); end
          %          size(examples)
          %          pause
          for v = 1:nFeatures
            if useCount & ~useSilence; if ~rem(v,1000); fprintf('%d ',v); end; end
            %           try
            neighbours = voxelsToNeighbours(v,1:numberOfNeighbours(v));

            %            neighbours
            %            examples(:,neighbours((end-1):end))'
            %            pause
            
            %           catch
            %             nFeatures
            %             size(numberOfNeighbours)
            %             size(voxelsToNeighbours)
            %             v
            %             numberOfNeighbours(v)
            %             voxelsToNeighbours(v,:)
            %            end
            %            [v,neighbours]
            tmp = examples(indicesTest,[v,neighbours]);
            tmp = tmp -  repmat(groupMean{lA}([v,neighbours]),nTest,1);
            tmp = tmp ./ repmat(groupStdv{lA}([v,neighbours]),nTest,1);
            tmp = -1 * (tmp .^2) / 2;
            tmp = sum(tmp - log(sqrt(2*pi)*repmat(groupStdv{lA}([v,neighbours]),nTest,1)),2);

            if usePriors; tmp = tmp + matrixOfPriors(:,lA); end

            % used to store scores directly, too much memory in a single structure
            %scores(indicesTest,v,lA) = sum(tmp,2);
            
            %scoresTest(:,v,lA) = sum(tmp,2); clear tmp;
            scoresTest(:,v,lA) = tmp; clear tmp;
          end
          if useCount & ~useSilence; fprintf('\n'); end

         case {'gnbtest'}
          % NOTE: groupStdv here is really variance
          if ~useSilence; fprintf('\tgroup %d: label %d: searchlight applying feature by feature (GNB test)\n',g,lA); end
          for v = 1:nFeatures
            if useCount & ~useSilence; if ~rem(v,1000); fprintf('%d ',v); end; end
            neighbours = voxelsToNeighbours(v,1:numberOfNeighbours(v));

            if 0
              tmp = examples(indicesTest,[v,neighbours]);
              tmp = (tmp -  repmat(groupMean{lA}([v,neighbours]),nTest,1)).^2;
              tmp = tmp ./ repmat(groupStdv{lA}([v,neighbours]),nTest,1);
              %tmp = sum(tmp - log(sqrt(2*pi)*repmat(groupStdv{lA}([v,neighbours]),nTest,1)),2);
              tmp = -0.5*sum(tmp,2);
            else
              tmp = 2*examples(indicesTest,[v,neighbours]).*repmat(groupMean{lA}([v,neighbours]),nTest,1) -...
                    repmat(groupMean{lA}([v,neighbours]),nTest,1) .^ 2 -...
                    examples(indicesTest,[v,neighbours]).^2;
              tmp = tmp ./ repmat(groupStdv{lA}([v,neighbours]),nTest,1);
              tmp = sum(tmp,2);
            end
            
            if usePriors; tmp = tmp + matrixOfPriors(:,lA); end

            scoresTest(:,v,lA) = tmp; clear tmp;
          end
          if useCount & ~useSilence; fprintf('\n'); end
          
        end
        
      end; % if searchlight
      
    end; % for over labels
    
    if useSearchlight
      if storeSearchlightCovarianceMatrices
        % no need for a loop over labels
        switch classifier
         case {'gaussianPooled','nbayesPooled','gnb','gnbtest','gnb_pooled'}
         otherwise
          fprintf('error: searchlight covariance matrix for non-pooled variance estimates does not make sense\n');return;
        end
        
        tmp = std(examples(indicesTrain,:),0,1);
        for l = 1:nLabels; groupStdv{l} = tmp; end

        for v = 1:nFeatures
          %if useCount & ~useSilence; if ~rem(v,1000); fprintf('%d ',v); end; end
          
          neighbours = voxelsToNeighbours(v,1:numberOfNeighbours(v));

          S = diag(groupStdv{1}([v,neighbours]));

          if storeSearchlightCovarianceMatrices
            % the covariance matrix
            optionalReturns.searchlightCovarianceMatrices(1:(numberOfNeighbours(v)+1),1:(numberOfNeighbours(v)+1),v) = optionalReturns.searchlightCovarianceMatrices(1:(numberOfNeighbours(v)+1),1:(numberOfNeighbours(v)+1),v) + S;          
          end
        end
      end
    end

%%%%% LDA/QDA, different loop     
    
   case {'lda','lda_shrinkage','lda_ridge','qda_shrinkage'}
    % this requires an entirely different loop

    %disp(examples(:,[(end-3) (end-2) (end-1) end])); pause
    
    % find class means and subtract them from respective training set points
    for l = 1:nLabels
      label = sortedLabelValues(l);
      
      indices{l} = find(labelsTrain == label); % training set examples with this label
      
      groupMean{l} = mean(examples(indicesTrain(indices{l}),:),1);

      examples(indicesTrain(indices{l}),:) = examples(indicesTrain(indices{l}),:)-repmat(groupMean{l},length(indices{l}),1);
    end; % for over labels
    
    %% main loop over features

    switch classifier
     case {'lda_shrinkage'}
      nb = size(voxelsToNeighbours,2)+1;
      W = zeros(nb,nb,nTrain);
     otherwise
    end

    switch classifier
     
     case {'lda','lda_ridge','lda_shrinkage'}
      for v = 1:nFeatures
        if useCount & ~useSilence; if ~rem(v,1000); fprintf('%d ',v); end; end

        neighbours = voxelsToNeighbours(v,1:numberOfNeighbours(v));

        % find the inverse of the covariance matrix (shared between classes)
        switch classifier
         case {'lda'}
          S = examples(indicesTrain,[v,neighbours])'*examples(indicesTrain,[v,neighbours])/(nTrain-1);
         case {'lda_ridge'}
          % lambda is fixed to a small number
          S = (examples(indicesTrain,[v,neighbours])'*examples(indicesTrain,[v,neighbours]) + eye(numberOfNeighbours(v)+1)*1e-10)/(nTrain-1);
         case {'lda_shrinkage'}
          S = cov_shrinkage(examples(indicesTrain,[v,neighbours]));
        end
        Sinv = inv(S);
        
        if storeSearchlightCovarianceMatrices
          optionalReturns.searchlightCovarianceMatrices(1:(numberOfNeighbours(v)+1),1:(numberOfNeighbours(v)+1),v) = optionalReturns.searchlightCovarianceMatrices(1:(numberOfNeighbours(v)+1),1:(numberOfNeighbours(v)+1),v) + S;
        
          %optionalReturns.searchlightMeans = zeros(nFeatures,27,nLabels);
          % the mean
          for l = 1:nLabels
            optionalReturns.searchlightMeans(v,1:(numberOfNeighbours(v)+1),l) = ...
                optionalReturns.searchlightMeans(v,1:(numberOfNeighbours(v)+1),l) +  groupMean{l}([v,neighbours]);
          end
        end
        
        % classify test examples
        for lA = 1:nLabels       
          tmp = examples(indicesTest,[v,neighbours]) - repmat(groupMean{lA}([v,neighbours]),nTest,1);
          scoresTest(:,v,lA) = diag(tmp*Sinv*tmp');
        end
      end; % loop over voxels

     case {'qda_shrinkage'}
      
      for v = 1:nFeatures
        if useCount & ~useSilence; if ~rem(v,1000); fprintf('%d ',v); end; end

        neighbours = voxelsToNeighbours(v,1:numberOfNeighbours(v));
        
        % classify test examples
        for lA = 1:nLabels       
          % find the <lA> class covariance matrix
          S = cov_shrinkage(examples(indicesTrainClass{lA},[v,neighbours]));
          Sinv = inv(S);
          tmp = examples(indicesTest,[v,neighbours]) - repmat(groupMean{lA}([v,neighbours]),nTest,1);
          % reverse sign
          scoresTest(:,v,lA) = - exp(-diag(tmp*Sinv*tmp')/2) / sqrt(det(S));
        end
      end; % loop over voxels
      
    end
    
    
    if useCount & ~useSilence; fprintf('\n'); end
    
    % replace class means
    for l = 1:nLabels
      examples(indicesTrain(indices{l}),:) = examples(indicesTrain(indices{l}),:)+repmat(groupMean{l},length(indices{l}),1);
    end; % for over labels
    
    scoresTest = exp(-0.5*scoresTest);
    
%%%%% all the SVM classifiers
    
   case {'svm_linear','svm_quadratic','svm_sigmoid','svm_rbf'}
    
    %
    % parameter processing (taken from the classifier wrappers):
    % - classifierLIBSVMwrapper
    % - classifierLogisticRegressionMulti
    %
    
    % defaults for options
    lambda = 1;
    gamma  = 1/nFeatures; % default value and term used in libsvm docs
    sigma  = sqrt(1/(2*gamma));
    coef   = 0;           % default value and term used in libsvm docs
    degree = 2; % comes in if using polynomial, not supported here yet
    crossvalidateLambda = 0;
    crossvalidateGamma  = 0;
    kernelDefined = 0;
    
    lcp = length(classifierParameters);
    idx = 1;
    while idx <= lcp
      argName = classifierParameters{idx}; idx = idx + 1;
      
      switch argName
       case {'kernel'}
        kernel = classifierParameters{idx}; idx = idx + 1;
        kernelDefined = 1;
       case {'lambda'}
        lambda = classifierParameters{idx}; idx = idx + 1;
        if isequal(lambda,'crossValidation')
          crossvalidateLambda = 1;
          if (isempty(groupLabels) | ~isrealmat(groupLabels) | (length(groupLabels)~=nExamples))
            fprintf('error: to select lambda with crossvalidation you need to provide group labels\n');pause;return;          
          end
        else
          if lambda <= 0; fprintf('error: lambda must be positive\n',this); pause;return; end
        end
        % parameter for kernels other than linear
       case {'degree'}
        if ~kernelDefined; fprintf('error: please specify kernel argument first\n');return;end
        degree = classifierParameters{idx}; idx = idx + 1;
       case {'gamma'}
        if ~kernelDefined; fprintf('error: please specify kernel argument first\n');return;end
        if isequal(kernel,'linear');
          fprintf('warning: no gamma with a linear kernel, gamma = 0\n'); 
          gamma = 0; idx = idx + 1;
        else
          gamma = classifierParameters{idx}; idx = idx + 1;
          if isequal(gamma,'crossValidation')
            crossvalidateGamma = 1;
            if (isempty(groupLabels) | ~isrealmat(groupLabels) | (length(groupLabels)~=nExamples))
              fprintf('error: to select gamma with crossvalidation you need to provide group labels\n');pause;return;          
            end
          else
            switch kernel
             case {'rbf'}
              sigma = sqrt(1/(2*gamma));
             otherwise
              % nothing else to derive
            end
          end
        end
       case {'sigma'}
        if ~kernelDefined; fprintf('error: please specify kernel argument first\n');return;end
        sigma = classifierParameters{idx}; idx = idx + 1;
        switch kernel
         case {'rbf'}
          gamma = 1/(2*sigma^2);
         otherwise
          fprintf('error: sigma argument only works in RBF kernels\n');return;
        end
       case {'coef'}
        if ~kernelDefined; fprintf('error: please specify kernel argument first\n');return;end
        coef = classifierParameters{idx}; idx = idx + 1;
       otherwise
        fprintf('error: unknown argument %s\n',argName); pause;return
      end; % switch on argname
    end; % loop over arguments

    % argument processing done inside classifierLIBSVMwrapper
    
    % defining a kernel argument string simplifies the rest of the code,
    % also define the parameter string for non-crossvalidated parameters
    
    parameterString = '';
    if ~crossvalidateLambda
      parameterString = sprintf('%s -c %s',parameterString,num2str(lambda));
    end
    if ~crossvalidateGamma & ~isequal(kernel,'linear')
      parameterString = sprintf('%s -g %s',parameterString,num2str(gamma));
    end
    
    % coefficient r is not cross-validated so might as well put it in kernelString
    switch kernel
     case {'linear'}
      kernelString = sprintf('-t 0 ');
     case {'quadratic','polynomial'}
      if isequal(kernel,'quadratic'); degree = 2; end
      kernelString = sprintf('-t 1 -d %d',degree);
      if coef > 0; kernelstring = sprintf('%s -r %s',kernelString,num2str(coef)); end
     case {'rbf'}
      kernelString = sprintf('-t 2 ');
      if coef > 0; kernelstring = sprintf('%s -r %s',kernelString,num2str(coef)); end
     case {'sigmoid'}
      kernelString = sprintf('-t 3 ');
      if coef > 0; kernelstring = sprintf('%s -r %s',kernelString,num2str(coef)); end
     otherwise
      fprintf('%s: unknown kernel type %s\n',this,kernel); pause;return;
    end

    if ~useSilence; fprintf('\tgroup %d: searchlight %s applying feature by feature\n',g,classifier); end
      
    %% if crossvalidation of either parameter was requested
    %% do a grid search (each combination of settings requires running a CV
    
    if crossvalidateLambda | crossvalidateGamma
      % if we need crossvalidation, set up the ranges
      fprintf('\tusing crossvalidation to find lambda/gamma\n');
      
      % wastes memory but cleans up the code, at least till it's debugged
      cvexamples    = examples(indicesTrain,:); ncvexamples = nTrain;
      cvlabels      = labels(indicesTrain);
      % use the original labels, in case group labels are just even/odd
      cvgroupLabels = groupLabelsOriginal(indicesTrain);
      %cvgroupLabels = groupLabels(indicesTrain);
      cvgroups = unique(cvgroupLabels); ncvgroups = length(cvgroups);
      
      if ncvgroups == 1
        fprintf('error: only 1 group labels within training set, consider using groupLabelsOriginal\n');pause;return;
      end
        
      bestcvresults = repmat(-1,1,nFeatures); % best accuracy thus far

      if crossvalidateLambda
        lambdaRange = [1 10 0.1 100 0.01];
      else
        lambdaRange = [lambda]; % fake range for outer loop
      end
      nl = length(lambdaRange);
      lambdas = repmat(lambdaRange(1),1,nFeatures);

      if crossvalidateGamma
        % range to try depends on the number of features
        gammaRange = [1/nFeatures, 10.^(ceil(log10(1/nFeatures)):log10(1))];
      else
        switch kernel
         case {'linear'}
          gammaRange = [0];
         otherwise
          gammaRange = [gamma]; % fake range for outer loop
        end
      end
      ng = length(gammaRange);
      gammas = repmat(gammaRange(1),1,nFeatures);
      
      % and grid search over all combinations of settings

      for ir = 1:nl
        lambda = lambdaRange(ir);

        baseString = sprintf(' -c %s',num2str(lambda));
        %fprintf('\ttesting lambda=%s\n',num2str(lambda));
        
        for ig = 1:ng
          gamma = gammaRange(ig);
          
          fprintf('\ttesting lambda=%s gamma=%s\n',num2str(lambda),num2str(gamma));
          
          if gamma
            parameterString = sprintf('%s -g %s',baseString,num2str(gamma));
          else
            % in linear kernel
            parameterString = baseString;
          end
          completeString  = sprintf('%s %s',kernelString,parameterString);
          
          %fprintf('\t\ttesting gamma=%s\n',num2str(gamma));
          
          cvresults = zeros(1,nFeatures);
          
          for icvg = 1:ncvgroups
            cvgroup = cvgroups(icvg);
            cvmask = (cvgroupLabels == cvgroup);
            
            %fprintf('\t\tcrossvalidation loop\n');
            cvindicesTest  = find( cvmask); cvnTest = length(cvindicesTest);
            cvindicesTrain = find(~cvmask); 

            zdummy = zeros(cvnTest,1); nPairs = nLabels*(nLabels-1)/2;
            tmp = zeros(cvnTest,nFeatures,nPairs); % hold predictions in this group

            % loop over voxels here
            for v = 1:nFeatures
              neighbours = voxelsToNeighbours(v,1:numberOfNeighbours(v));

              % now run the loop
              cvModel = svmtrain(cvlabels(cvindicesTrain),cvexamples(cvindicesTrain,[v,neighbours]),completeString);              
              
              [discard1,discard2,tmp(:,v,:)] = svmpredict(zdummy,cvexamples(cvindicesTest,[v,neighbours]),cvModel);
            end

            % transform decision values into vote counts for each class
            votes = zeros(cvnTest,nFeatures,nLabels);
            ip = 1;
            for ic1 = 1:(nLabels-1)
              for ic2 = (ic1+1):nLabels
                vmask = (tmp(:,:,ip) > 0);

                % positive are votes for the first class in the pair
                votes(:,:,ic1) = votes(:,:,ic1) +  vmask;
                % negative are votes for the second class in the pair
                votes(:,:,ic2) = votes(:,:,ic2) + ~vmask;
                % ties divide the vote
                vmask = (tmp(:,:,ip) == 0)/2;
                votes(:,:,ic1) = votes(:,:,ic1) + vmask;
                votes(:,:,ic2) = votes(:,:,ic2) + vmask;

                ip = ip + 1;
              end
            end
            
            [discard,tmp] = max(votes,[],3);
            
            cvresults = cvresults + sum(tmp == repmat(cvlabels(cvindicesTest),1,nFeatures),1);
            clear votes tmp discard1 discard2 discard;
            
          end; % loop over groups
          cvresults = cvresults / ncvexamples; % accuracy
          
          % update lambda/gamma for voxels where current values gave better results
          indicesChanged = find(cvresults > bestcvresults);
          bestcvresults(indicesChanged) = cvresults(indicesChanged);
          lambdas(indicesChanged) = lambda;
          gammas(indicesChanged)  = gamma;
          
        end; % loop over gamma
      end; % loop over lambda
      
    else
      % no need to crossvalidate to find lambda/gamma
      lambdas = repmat(lambda,1,nFeatures); gammas = repmat(gamma,1,nFeatures);
      fprintf('\tusing given lambda/gamma (%s/%s)\n',num2str(lambda),num2str(gamma));
    end; % if cross-validating
    
    % store lambda/gamma used in each voxel
    optionalReturns.lambdas = lambdas;
    optionalReturns.gammas  = gammas;
    
    %% now train/test classifier on the full training data
    %% using estimated or provided lambdas/gammas
    
    fprintf('\tlearning full training set SVM\n');
    zdummy = zeros(nTest,1); nPairs = nLabels*(nLabels-1)/2;

    % hold decision values in this group
    %tmp = zeros(nTest,nFeatures,nLabels);
    tmp = zeros(nTest,nFeatures,nPairs); % works for two classes
    
    %tmplabels = zeros(nTest,nFeatures); % use for multi-label debugging below

    caralho = find(lambdas <= 0);
    if ~isempty(caralho)
      fprintf('error! there are lambdas <= 0\n'); 
      lambdas(caralho)
      return;
    end

    % precompute argument strings to make this cleaner
    switch kernel
     case {'linear'}
      for v = 1:nFeatures; stringvoxel{v} = sprintf('%s -c %s',kernelString,num2str(lambdas(v)));end
     otherwise
      for v = 1:nFeatures; stringvoxel{v} = sprintf('%s -c %s -g %s',kernelString,num2str(lambdas(v)),num2str(gammas(v)));end
    end
    
    for v = 1:nFeatures
      neighbours = voxelsToNeighbours(v,1:numberOfNeighbours(v));
      if useCount & ~useSilence; if ~rem(v,1000); fprintf('%d ',v); end; end        

      % train
      model = svmtrain(labels(indicesTrain),examples(indicesTrain,[v,neighbours]),stringvoxel{v});
      
      [discard1,discard2,tmp(:,v,:)] = svmpredict(zdummy,examples(indicesTest,[v,neighbours]),model);
      %[tmplabels(:,v),discard2,tmp(:,v,:)] = svmpredict(zdummy,examples(indicesTest,[v,neighbours]),model);
    end; % loop over voxels

    % what to do with decision values depends on the # labels
    % WARNING: assumes the label of class 1 is less than the label of class 2
    % (see below)
    
    % transform decision values into vote counts for each class
    votes = zeros(nTest,nFeatures,nLabels);
    ip = 1;
    for ic1 = 1:(nLabels-1)
      for ic2 = (ic1+1):nLabels
        vmask = (tmp(:,:,ip) > 0);

        % positive are votes for the first class in the pair
        votes(:,:,ic1) = votes(:,:,ic1) +  vmask;
        % negative are votes for the second class in the pair
        votes(:,:,ic2) = votes(:,:,ic2) + ~vmask;
        % ties divide the vote
        vmask = (tmp(:,:,ip) == 0)/2;
        votes(:,:,ic1) = votes(:,:,ic1) + vmask;
        votes(:,:,ic2) = votes(:,:,ic2) + vmask;

        ip = ip + 1;
      end
    end
    scoresTest = votes; clear votes tmp discard1 discard2;
    
    % check that we are replicating the voting mechanism correctly
    % make sure you uncomment tmplabels definition above,
    % then the line with tmplabels(:,v)
    % and then uncomment this:
    %[discard,maxpos] = max(scoresTest,[],3);
    %isequal(tmplabels,sortedLabelValues(maxpos))
    %pause
    
    %    if model.Label(1) < model.Label(2)
    %      % labels are in the right order
    %      scoresTest(:,:,1) = tmp;
    %      scoresTest(:,:,2) = -tmp;
    %    else
    %      % labels are flipped
    %      scoresTest(:,:,1) = -tmp;
    %      scoresTest(:,:,2) = tmp;
    %    end
    fprintf('\n');
    
%%%% nearest neighbour classifiers
    
   case {'knn'}

    % process arguments, as the algorithms are different for
    % k=0 - distance to class means
    % k=1 - nearest neighbour
    % k>1 - k-nn proper
    
    % defaults and processing of overrides
    k = 0;
    similarityMeasure = 'correlation';
    
    lcp = length(classifierParameters);
    idx = 1;
    while idx <= lcp
      argval = classifierParameters{idx}; idx = idx + 1;
      
      switch argval
       case {'k'}
        k = classifierParameters{idx}; idx = idx + 1;
        k = ceil(k); if k < 0; fprintf('error: knn: k < 0\n');return;end
        
       case {'similarityMeasure'}
        similarityMeasure = classifierParameters{idx}; idx = idx + 1;
       otherwise
        fprintf('error: knn: unknown argument %s\n',argval);return;
      end
    end

    %% main loop is over classes
    
    % ancillary matrices
    switch k
     case {0}
      tmpTrain = exampleMeanTrainClass';
      tmpTest  = examples(indicesTest,:)';
     case {1}
      tmpTrain = examples(indicesTrain,:)';
      tmpTest  = examples(indicesTest,:)';
    end
    
    switch similarityMeasure
     case {'correlation'}
      measurenum = 0;
     case {'euclidean'}
      measurenum = 1;
    end
    
    tmp = (1:nTest)';
    
    % compute similarity measure for each neighbourhood
    for v = 1:nFeatures
      neighbours = voxelsToNeighbours(v,1:numberOfNeighbours(v));
      if useCount & ~useSilence; if ~rem(v,1000); fprintf('%d ',v); end; end        

      [St,indicesMostSimilar,valuesMostSimilar] = rowiseSimilarity(tmpTest([v,neighbours],:),tmpTrain([v,neighbours],:),measurenum);
      
      S = St'; % #test x #classes (k=0) or x #train (k=1)
      % indicesMostSimilar - #test x 1 - entries are in 1:#classes range (k=0) or 1:#train (k=1)
      
      % all scores are 0 except 1 for the predicted class
      template = zeros(nTest,nLabels);
      switch k
       case {0}
        indices = sub2ind([nTest,nLabels],tmp,indicesMostSimilar);

%        if numberOfNeighbours(v) == 0
%          size(indicesMostSimilar)
%          nLabels
%          [indicesMostSimilar,valuesMostSimilar]
%          neighbours
%          pause
%        end
          
        template(indices) = 1;
        scoresTest(:,v,:) = template;
       case {1}
        % labels are already 1:#classes
        indices = sub2ind([nTest,nLabels],tmp,labels(indicesTrain(indicesMostSimilar)));
        template(indices) = 1;
        scoresTest(:,v,:) = template;
      end

      %disp([labels(indicesTest),template]);pause
      
    end; % for over features
      
%      switch classifier
%       case {'lr_L1'}
%        model = classifierLogisticRegressionMulti(examples(indicesTrain,[v,neighbours]),labels(indicesTrain),'regularization',{'L1',lambda},'useSilence');
%       case {'lr_L2'}
%        model = classifierLogisticRegressionMulti(examples(indicesTrain,[v,neighbours]),labels(indicesTrain),'regularization',{'L2',lambda},'useSilence');
%      end
%      decisions = [ones(nTest,1),examples(indicesTest,[v,neighbours])] * model.W;
%      %scoresTest(:,v,:) = decisions;
%      scoresTest(:,v,1) = decisions(:,1);
%      scoresTest(:,v,2) = decisions(:,2);    
%    end; % loop over voxels
    
    
    
%%%% general purpose classifiers 
    
   case {'lr_L2','lr_L1'}

    lambda = 1;
      
    for v = 1:nFeatures
      neighbours = voxelsToNeighbours(v,1:numberOfNeighbours(v));
      if useCount & ~useSilence; if ~rem(v,1000); fprintf('%d ',v); end; end        
      
      switch classifier
       case {'lr_L1'}
        model = classifierLogisticRegressionMulti(examples(indicesTrain,[v,neighbours]),labels(indicesTrain),'regularization',{'L1',lambda},'useSilence');
       case {'lr_L2'}
        model = classifierLogisticRegressionMulti(examples(indicesTrain,[v,neighbours]),labels(indicesTrain),'regularization',{'L2',lambda},'useSilence');
      end
      decisions = [ones(nTest,1),examples(indicesTest,[v,neighbours])] * model.W;
      %scoresTest(:,v,:) = decisions;
      scoresTest(:,v,1) = decisions(:,1);
      scoresTest(:,v,2) = decisions(:,2);    
    end; % loop over voxels      

   otherwise
    fprintf('error: should not get here with %s\n',classifier);pause;return;
    
  end; % switch over classifiers
    
  %
  % for each feature, compute the error measure (error or average rank) in each example
  %
  
  template = repmat(labelsTest,1,nFeatures);
  
  switch measure
   case {'error'}
    clear labelRankings;
    
    [maxScores,maxIndices] = max(scoresTest,[],3);  clear scoresTest;
    predictedLabels = sortedLabelValues( maxIndices );
    
    measurePerExample(indicesTest,:) = (predictedLabels ~= template);
    
    
   case {'averageRank'}
    % last slice of sortedIndices contains z-coordinates of highest scores
    
    [sortedScores,sortedIndices] = sort(scoresTest,3,'descend');  clear scoresTest;
    labelRankings   = sortedLabelValues( sortedIndices ); clear sortedIndices;
    predictedLabels = labelRankings(:,:,1);
    
    % labelRankings has, for each example (x) and feature (y)
    % a ranking (z) of the labels by score (highest is first)
    
    % Same dimensions, 0 everywhere except 1 where the true labels are
    % in the ranking. The position of the 1 in each ranking is the position
    % of the correct label in that ranking
    for l = 1:nLabels
      labelRankings(:,:,l) = (labelRankings(:,:,l) == template);
    end
    %labelRankingsMask = (labelRankings == label); clear labelRankings;
      
    % finds the position of the 1
    [dummy,correctLabelRank] = max(labelRankings,[],3);
    clear labelRankings;
    
    % compute normalized rank (varies between 0 (perfect) and 1)
    measurePerExample(indicesTest,:) = (correctLabelRank - 1) / (nLabels - 1);
  end
  
  %% if required, update the local confusion matrix
  
  if storeLocalConfusionMatrices
    ncc = nLabels*nLabels;
    indexMatrix = repmat((0:ncc:(ncc*(nFeatures-1))),nTest,1);
    indexMatrix = indexMatrix + (predictedLabels-1)*nLabels + template;
    
    for i = 1:(nFeatures*nTest)
      optionalReturns.localConfusionMatrices(indexMatrix(i)) = optionalReturns.localConfusionMatrices(indexMatrix(i)) + 1;
    end
  
  end
  
end; % for over groups

if storeSearchlightCovarianceMatrices
  % up till now we have been adding these across groups
  % now divide by the number of groups to average
  optionalReturns.searchlightCovarianceMatrices = optionalReturns.searchlightCovarianceMatrices /nGroups;
  optionalReturns.searchlightMeans = optionalReturns.searchlightMeans / nGroups;
end

%
% compute errors
%

measureOverall  = zeros(1,nFeatures);
measurePerClass = cell(nLabels,1);

for lB = 1:nLabels
  %% for examples of class <lB>, sort label scores
  label = sortedLabelValues(lB);
  %    measurePerClass{lB} = sum(measurePerExample(indicesLabel{lB},:),1)/nPerLabel;
  measurePerClass{lB} = sum(measurePerExample(indicesLabel{lB},:),1)/nPerLabel(lB);
end

% average measure over all examples
measureOverall = sum(measurePerExample,1) / nExamples;

% make a binary matrix (#examples x #voxels) with correct/incorrect
if storeCorrectMatrix
  optionalReturns.correctMatrix = measurePerExample;
end

if storeAverageClassMeansStdvs
  for l = 1:nLabels
    optionalReturns.averageClassMeans(l,:) = optionalReturns.averageClassMeans(l,:) / nGroups;
    optionalReturns.averageClassStdvs(l,:) = optionalReturns.averageClassStdvs(l,:) / nGroups;
  end
end

% for debugging
%  tmpAccuracy = zeros(1,nFeatures);
%  for v = 1:nFeatures
%    tmpAccuracy(v) = 1-(sum(diag(optionalReturns.localConfusionMatrices(:,:,v)))/nExamples);
%  end
%  tmp = reshape(optionalReturns.localConfusionMatrices,[nLabels*nLabels nFeatures]);


  
%% finally, sort the features
[sortedMeasure,sortedFeatures] = sort(measureOverall);





%
% Takes a matrix and computes mean across all rows except one.
% The computation is done for every row in turn (see line below
% for details)
%
% In:
% - <m> x <n> matrix to perform the computation on
%
% Out:
% - <m> x <n> matrix, where row <r> contains the mean over all rows
%   except row <r>
%
 
function [matrixOfMeans] = aux_meansMinusOneRow( m )
 
[nrows,ncols] = size(m);
matrixOfMeans = zeros(nrows,ncols);
 
matrixSum = sum(m,1);
 
for r = 1:nrows
  matrixOfMeans(r,:) = matrixSum - m(r,:);
end
 
matrixOfMeans = matrixOfMeans ./ (nrows-1);


%
% Takes a matrix and computes std across all rows except one.
% The computation is done for every row in turn (see line below
% for details)
%
% In:
% - <m> x <n> matrix to perform the computation on
%
% Out:
% - <m> x <n> matrix, where row <r> contains the standard deviation
%   over all rows, except row <r>
%
                                                                                                 
function [matrixOfStdevs] = aux_stdevsMinusOneRow( m )
                                                                                                 
[nrows,ncols] = size(m);
matrixOfStdevs = zeros(nrows,ncols);
                                                                                                 
m1 = m;
m2 = m .^ 2;
                                                                                                 
matrixOfMeans1 = aux_meansMinusOneRow(m1);
matrixOfMeans2 = aux_meansMinusOneRow(m2);
                                                                                                 
matrixOfStdevs = matrixOfMeans2 - (matrixOfMeans1 .^ 2);
matrixOfStdevs = sqrt(matrixOfStdevs);




