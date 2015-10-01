%
% Takes several pvalue maps and
% - finds which voxels are significantly different from chance in each map
%   (using FDR)
% - across voxels significant in any map, compares all pairs of maps
%
% In:
% - accuracyMaps #1 - #maps x #voxels
% - pvalueMaps   #1 - #maps x #voxels
% - optional
%   - accuracyMaps #2 - #maps x #voxels
%   - pvalueMaps   #1 - #maps x #voxels
% - optional (have to come after all the maps)
%   - 'q' - FDR q-value for each map (defaults to 0.01)
%   - 'p' - p-value for the post-hoc test (defaults to 0.01)
%
% Out:
% - significantPairs   - #maps x #maps binary array with which pairs were significantly different
% - significantPvalues - same, but now each cell has the p-value multiplied by the sign
%                        (for plotting purposes)
%

function [results] = compareInformationMaps(varargin)

%% process arguments

if nargin < 2;
  fprintf('syntax: compareInformationMaps(<accuracyMaps>,<pvalueMaps>,[see help]\n');return;
end
  
accuracyMapsA = varargin{1}; [nmapsA,nvoxels] = size(accuracyMapsA);
pvalueMapsA   = varargin{2};
if ~isequal(size(pvalueMapsA),[nmapsA,nvoxels])
  fprintf('error: accuracyMapsA and pvalueMapsA must be the same size\n');return;
end

idx = 3;
if nargin > 2
  if isrealmat(varargin{idx})
    % comparison is with a different set of accuracy maps (will crash if not two args)
    accuracyMapsB = varargin{idx}; idx = idx + 1;
    pvalueMapsB   = varargin{idx}; idx = idx + 1;
    sameDataset = 0;
  else
    % same sets, will create a reference rather than duplicate the matrix
    accuracyMapsB = accuracyMapsA;
    pvalueMapsB   = pvalueMapsA;
    sameDataset = 1;
  end
end
[nmapsB,nvoxelsB] = size(accuracyMapsB);
if ~isequal(size(pvalueMapsB),[nmapsB,nvoxelsB])
  fprintf('error: accuracyMapsB and pvalueMapsB must be the same size\n');return;
end
if nvoxels ~= nvoxelsB
  fprintf('error: accuracyMapsA and accuracyMapsB must have the same number of voxels\n');return;
end

comparisonToMake     = 'spatialCorrelation';
multipleComparison   = 'fdr';
comparisonTest       = 'signedrank';
comparisonCorrection = 'hsd';

while idx <= nargin
  argval = varargin{idx}; idx = idx + 1;
    
  switch argval
   case {'comparisonToMake'}
    comparisonToMake = varargin{idx}; idx = idx + 1;
    
    switch comparisonToMake
     case {'spatialCorrelation'}
      % no need for more arguments
     case {'significanceMaps'}
      % requires two more, test and significance level
      multipleComparison = varargin{idx}; idx = idx + 1;
      significanceLevel  = varargin{idx}; idx = idx + 1;
     
     case {'significantAccuracy'}
      % requires two more, p for post-hoc test and q for significant test
      switch varargin{idx}
       case {'default'}
        multipleComparison = 'fdr';
        significanceLevel    = 0.01;
        comparisonTest       = 'signedrank';
        comparisonCorrection = 'hsd';
        comparisonLevel      = 0.01;
        idx = idx + 1;
       otherwise
        % assumes all of them will be specified
        if (idx + 4) > nargin
          fprintf('error: significantAccuracy requires four arguments: <multipleComparison=*fdr*|bonferroni>,<significance level>,<comparisonTest=*signedrank*|ttest>,<comparison level>\n'); 
        else
          multipleComparison   = varargin{idx}; idx = idx + 1;
          significanceLevel    = varargin{idx}; idx = idx + 1;
          comparisonTest       = varargin{idx}; idx = idx + 1;
          comparisonCorrection = varargin{idx}; idx = idx + 1;
          comparisonLevel      = varargin{idx}; idx = idx + 1;
        end
      end
    end; % of switch on compariso
  
   otherwise
    fprintf('error: unknown argument %s\n',argval);return;
    
  end; % of switch on argval
end

switch multipleComparison
 case {'fdr'}
 otherwise
  fprintf('error: %s multiple comparison is not supported\n',multipleComparison);return;
end
  
switch comparisonTest
 case {'signedrank'}
 case {'ttest'}
 otherwise
  fprintf('error: %s accuracy comparison is not supported\n',comparisonTest);return;  
end

%
% Compute significance in each map
%

S = zeros(nmapsA,nmapsB);

switch comparisonToMake

  %
  % Spatial correlation
  %
 case {'spatialCorrelation'}
  
  if sameDataset
    results.S = computeExampleSimilarity(accuracyMapsA,'correlation');
  else
    results.S = computeExampleSimilarity(accuracyMapsA,accuracyMapsB,'correlation');
  end
    
  %
  % Significant accuracy or just significance maps
  %
 case {'significanceMaps','significantAccuracy'}

  % determine significant voxels in each map of the two sets
  switch multipleComparison
   case {'fdr'}
    [thresholdID,thresholdN] = computeFDR(pvalueMapsA',significanceLevel);
    %threshold = thresholdID';
    threshold = thresholdN';
    significanceMapsA = (pvalueMapsA <= repmat(threshold,1,nvoxels));
    
    if sameDataset
      significanceMapsB = significanceMapsA;
    else
      [thresholdID,thresholdN] = computeFDR(pvalueMapsB',significanceLevel);
      %threshold = thresholdID';
      threshold = thresholdN';
      significanceMapsB = (pvalueMapsB <= repmat(threshold,1,nvoxels));
    end
  end

  % find voxels that are significant on any map
  significanceMaskA = max(significanceMapsA,[],1);
  significanceIndicesA = find(significanceMaskA);
  
  if sameDataset
    significanceMask = significanceMaskA;
  else
    significanceMaskB = max(significanceMapsB,[],1);
    significanceIndicesB = find(significanceMaskB);
    significanceMask = significanceMaskA | significanceMaskB;
    significanceIndices = find(significanceMask);
  end
  significanceIndices = find(significanceMask); nsig = length(significanceIndices);
    
  switch comparisonToMake
   case {'significanceMaps'}
    % all that we need to return are the significance maps
    if sameDataset
      results.significanceMapsA = significanceMapsA;
    else
      results.significanceMapsA = significanceMapsA;
      results.significanceMapsB = significanceMapsB;      
    end
    results.significanceIndices = significanceIndices;
    return;
   otherwise
    % keep going
  end

  % for each pair of maps

  npairs = nmapsA*nmapsB;

  if 1
    % needs data shaped as #voxels x #classifiers to compare
    accuracySignificant = accuracyMapsA(:,significanceIndices)';
    
    [differenceMatrix,significanceMatrix,rankMatrix,differencePair] = compareResultsFromMultipleClassifiers(accuracySignificant);
    results.differenceMatrix   = differenceMatrix;
    results.significanceMatrix = significanceMatrix;
    results.rankMatrix         = rankMatrix;
  else
    differencePair = zeros(npairs,nsig); % differences for each pair over global significant
    D = zeros(nmapsA,nmapsB);            % median paired difference
    P = zeros(nmapsA,nmapsB);            % p-value
    
    ip = 1;
    for ia = 1:nmapsA
      mapA = accuracyMapsA(ia,significanceIndices);
      for ib = 1:nmapsB
        mapB = accuracyMapsB(ib,significanceIndices);
        
        % apply the test
        differences = (mapA - mapB)';
        switch comparisonTest
         case {'signedrank'}
          [pvalue,h] = signrank(mapA',mapB');
          P(ia,ib) = pvalue;
         case {'ttest'}
          [h,pvalue] = ttest(mapA',mapB');
          P(ia,ib) = pvalue;
        end
        differencePair(ip,:) = differences;
        D(ia,ib) = median(differences);

        if sameDataset & (ia == ib); P(ia,ib) = 0; end
        
        ip = ip + 1;
      end
    end

    tmp = triu(P);
    [thresholdID,thresholdN] = fdr(tmp(find(tmp>0)),comparisonLevel);
%    sum(tmp(find(tmp>0))<=thresholdID)
%    pause
    % pack returns
    %results.significanceMatrix = P <= comparisonLevel;   % uncorrected
    results.significanceMatrix = P <= thresholdN;
    results.differenceMatrix = D;
    results.P = P;
  end
    
  results.significanceMapsA   = significanceMapsA;
  results.significanceMapsB   = significanceMapsB;
  results.significanceIndices = significanceIndices;
  results.differencePair      = differencePair;
  
  % for each pair of maps, run t-test or signed rank test on accuracy
  % (check for paired vs not)  
end


function [] = test()

nmaps = 5; nvoxels = 100;

accuracyMaps = rand(nmaps,nvoxels);
pvalueMaps   = rand(nmaps,nvoxels);

[results] = compareInformationMaps(accuracyMaps,pvalueMaps,'comparisonToMake','spatialCorrelation');
