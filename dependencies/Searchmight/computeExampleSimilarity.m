% Computes distance/similarity measures between the examples in one set or examples in two sets
%
% Input:
%
% - examples - #examples x #voxels (or features)
% - measure  - similarity/distance measure, can be: correlation | cosine | euclidean
%
% or
%
% - examplesA
% - examplesB
% - measure
%
% Output:
%
% - measure - #examples x #examples
%
% or
%
% - measure - #examplesA x #examplesB
%
% Notes:
% - measures:
%   - euclidean
%   - L1
%   - correlation
%   - cosine
%   - binary
%     - binary_incommon -
%     - binary_allones  -
%
%
% Examples:
%
% History:
% 2009 Feb 6 - fpereira@princeton.edu - created from earlier code
%

function [measureMatrix] = compute_exampleSimilarity(varargin)

%
% process arguments and data
%

this = 'compute_exampleSimilarity';

if nargin >= 2
  if nargin == 2
    examplesA  = varargin{1};
    measure    = varargin{2};
    examplesB  = examplesA; % this is just a pointer as long as examplesB is unchanged
    sameExamples = 1;
  else
    examplesA  = varargin{1};
    examplesB  = varargin{2};
    measure    = varargin{3};
    sameExamples = 0;
  end

  % possible extra arguments

  findClosestK = 0; % by default, distances to all the neighbours
  if nargin > 3
    idx = 4;
    while idx <= nargin
      argval = varargin{idx}; idx = idx + 1;
      switch argval
       case {'findClosestK'}
        findClosestK = varargin{idx}; idx = idx + 1;
       otherwise
        fprintf('%s: error: unknown argument %s\n',this,argval); return;
      end
    end
  end
  
  [nExamplesA,nVoxelsA] = size(examplesA); nA = nExamplesA;
  [nExamplesB,nVoxelsB] = size(examplesB); nB = nExamplesB;
  
  if (nVoxelsA ~= nVoxelsB);
    fprintf('%s: error: matrices need to have the same number of columns (features)\n',this);return;
  else
    nVoxels = nVoxelsA; m = nVoxels;
  end
  
  switch measure
   case {'euclidean'}
   case {'L1'}
   case {'cosine'}
   case {'correlation'}
   case {'binary_incommon'}
   case {'binary_allones'}
   case {'hellinger'}
   case {'KLdivergence'}
   case {'adhoc_percentsAndCounts'}
   otherwise
    fprintf('%s: error: measure %s is not supported\n',this,measure);return;
  end
else
  fprintf('syntax:\n');return
end

%
% Compute distances/similarities
%

%fprintf('%s: computing distances\n',this);

%% examples can be normalized ahead of time to make further calculation simpler

switch measure
 case {'KLdivergence'}
  fprintf('KLdivergence(P|Q): assuming A is P and B is Q\n');
  %fprintf('STOP! This has not been tested yet\n');pause;return;
  
  %  disp([sum(examplesA,2),sum(examplesB,2)]);pause
  sumA = sum(sum(examplesA,2)-ones(nA,1));
  sumB = sum(sum(examplesB,2)-ones(nB,1));
  if (sumA > 0.001) | (sumB > 0.001)
    fprintf('error: to use KL divergence all examples need to sum to 1\n');return;
  end
  
  % a lot of entries are 0, so make log(them) 0 rather than NaN
  indicesZero = find(examplesA==0);
  lexamplesA = log(examplesA); lexamplesA(indicesZero) = 0;
  indicesZero = find(examplesB==0);
  lexamplesB = log(examplesB); lexamplesB(indicesZero) = 0;

 case {'euclidean','L1','binary_incommon'}
  % nothing to do
 case {'correlation'}
  % make each image mean 0 and standard deviation one
  examplesA = examplesA -  repmat(mean(examplesA,2),1,m);
  examplesA = examplesA ./ repmat(std(examplesA,0,2),1,m);
  if sameExamples
    clear examplesB; examplesB = examplesA;
  else
    examplesB = examplesB -  repmat(mean(examplesB,2),1,m);
    examplesB = examplesB ./ repmat(std(examplesB,0,2),1,m);
  end

 case {'cosine'}
  % make each image L2 norm 1
  examplesA = examplesA ./ repmat(sqrt(sum(examplesA.^2,2)),1,m);
  if sameExamples
    clear examplesB; examplesB = examplesA;  
  else
    examplesB = examplesB ./ repmat(sqrt(sum(examplesB.^2,2)),1,m);
  end
  
 case {'hellinger'}
  % take the square root of all the examples coming in, as the formula does that
  examplesA = sqrt(examplesA);
  examplesB = sqrt(examplesB);

 case {'adhoc_percentsAndCounts'}
  % examplesA are percentages, examplesB are counts
  % each exampleA gets turned into fake counts, by multiplying by each exampleB total
  examplesBtotals = sum(examplesB,2);
examplesBorg = examplesB;
  % make each image mean 0 and standard deviation one for examplesB
  examplesB = examplesB -  repmat(mean(examplesB,2),1,m);
  examplesB = examplesB ./ repmat(std(examplesB,0,2),1,m);
end

%% computation proper

if findClosestK == 0
  % calculate distance/similarity from each example in A to all those in B
  measureMatrix = zeros(nA,nB);
else
  % if nA or nB are very large, measureMatrix will be enormous, this lets us
  % store just the few closest/furtherst distances
  measureMatrix = sparse(nA,nB);
end

for e = 1:nA
  switch measure
   case {'euclidean'}
    tmp = sum((repmat(examplesA(e,:),nB,1)-examplesB).^2,2);
   case {'KLdivergence'}
    tmp = (repmat(lexamplesA(e,:),nB,1) - lexamplesB) .* repmat(examplesA(e,:),nB,1);
    tmp = sum(tmp,2);
   case {'L1'}
    tmp = sum(abs(repmat(examplesA(e,:),nB,1)-examplesB),2);
   case {'correlation'}
    tmp = sum(repmat(examplesA(e,:),nB,1).*examplesB,2)/(m-1);
   case {'cosine'}
    tmp = sum(repmat(examplesA(e,:),nB,1).*examplesB,2);
   case {'binary_incommon'}
    % considers only the features that either example in a pair have
    % (note that this means that the features considered might be different
    %  when comparing A to B and A to C)
    tmp = repmat(examplesA(e,:),nB,1);
    maskOR  = tmp | examplesB; sumOR  = sum(maskOR,2);  clear maskOR;
    maskAND = tmp & examplesB; sumAND = sum(maskAND,2); clear maskAND;
    tmp = sumAND ./ sumOR;
   case {'binary_allones'}
    featuresToConsider = find(sum(examplesA,1) + sum(examplesB,1)); nf = length(featuresToConsider)
    tmp = repmat(examplesA(e,:),nB,1);
    
    maskAND = tmp(:,featuresToConsider) & examplesB(:,featuresToConsider); sumAND = sum(maskAND,2); clear maskAND;
    tmp = sumAND / nf;
    % most values will be high, because of all the 0s in common, normalize later
   case {'hellinger'}
    tmp = sum((repmat(examplesA(e,:),nB,1)-examplesB).^2,2)/2;
  
   case {'adhoc_percentsAndCounts'}
    % create fake counts and normalize them for correlation
    tmpc = round(repmat(examplesA(e,:),nB,1) .* repmat(examplesBtotals,1,m));
    tmpc = tmpc -  repmat(mean(tmpc,2),1,m);
    tmpc = tmpc ./ repmat(std(tmpc,0,2),1,m);
    tmp = sum(tmpc.*examplesB,2)/(m-1);
  end

  % depending on what is required, store straight or find few highest/sallest values
  if findClosestK == 0
    % store distance of this example to all the examples in matrix B
    measureMatrix(e,:) = tmp;
  else
    % find the closest/highest examples in B
    switch measure
     case {'euclidean','L1','binary_incommon','binary_allones','hellinger','KLdivergence'}
      [betterValues,betterFeatures] = sort(tmp);
     case {'correlation','adhoc_percentsAndCounts'}
      [betterValues,betterFeatures] = sort(tmp,'descend');
    end
    % make 0 values different from 0, the background in the matrix, by adding a very small number
    %measureMatrix(e,betterFeatures(1:findClosestK)) = betterValues(1:findClosestK) + 1e-8;
    measureMatrix(e,betterFeatures(1:findClosestK)) = betterValues(1:findClosestK);
  end
  
end

%% some post-facto normalization depending on the measures computed for all examples

switch measure
 case {'binary_allones'}
  % make it [0,1], with 0 being the lowest score and 1 being the highest
  minval = min(measureMatrix(:)); maxval = max(measureMatrix(:));
  measureMatrix = (measureMatrix - minval) ./ (maxval - minval);
  
 otherwise
  % no need
end
