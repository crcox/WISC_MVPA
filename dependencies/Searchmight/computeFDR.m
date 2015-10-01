% Computes the appropriate threshold to use over the pvalues
% resulting from multiple tests of a null hypothesis,
% in order to have a particular False Discovery Rate (FDR).
%
% This is the proportion of false positives (incorrect rejections
% of the null hypothesis) among those tests for which the null
% hypothesis is rejected. In voxel terms, this is the expected proportion
% of voxels declared active (using the threshold) that is not
% active. 
%
% Note that "expected" here means that, were one to replicate the
% experiment many times, the average FDR over those replications
% would be no bigger than the value for which the threshold was
% determined. For ANY PARTICULAR data analysis, the actual FDR
% might be LARGER than that value.
%
% For more details, see "Controlling FDR in Neuroimaging", by
% Genovese, Lazar and Nichols (in NeuroImage).
%
% Nichols' FDR web page says, apropos of this:
%
% Note on the preprint: The preprint reflected the current results at the time of submission. Since then, the independence FDR result of 1995 (BH 1995; in the above paper, "c(V)=1" on pg 9) has been shown to be more general (BY 2001). In particular, a technical condition known as positive regression dependency is sufficient to use the 1995 result. Examples of this condition include multivariate normal data where the covariance matrix has all positive elements; this seems a reasonable assumption for imaging, and hence the code above uses the c(V)=1 result. 
%
% so one should use the thresholdID output (which uses c(v)=1)
%
% In:
% - a vector of pvalues or a matrix whose columns are vectors of pvalues
% - a false discovery rate (e.g. 0.01, 0.05)
%
% Out:
% - the two thresholds described in the paper
%   - thresholdID - threshold based on independence or positive dependence
%   - thresholdN  - nonparametric pvalue threshold
%
% Dependencies:
%
% History
% - 2009 November 4 - created - fpereira@princeton.edu
%
%

function [thresholdID,thresholdN] = computeFDR( varargin )
  
  l = length(varargin);
  if l < 1
    fprintf('syntax: computeFDR(<pvalue vector or matrix with pvalue vector columns>,<fdr>\n');
    fprintf('please see code comment header for more info\n');
    return;
  end

  pValues = varargin{1}; [n,ncols] = size(pValues);
  fdr     = varargin{2};

  if n == 1; pValues = pValues'; n = ncols; ncols = 1; end
  
  significant = zeros(n,1);

  [pValuesSorted,sortedIndices] = sort(pValues,1);

  constantID = 1;
  constantN  = sum(1./(1:n));
  
  ratioID = fdr / (n*constantID);
  ratioN  = fdr / (n*constantN);

  pos = repmat((1:n)',1,ncols);
  maskID = (pValuesSorted ./ pos) <= ratioID;
  [maxrowID] = max(maskID .* pos,[],1);
  maskN = (pValuesSorted ./ pos) <= ratioN;
  [maxrowN] = max(maskN .* pos,[],1);

  thresholdID = zeros(1,ncols);
  thresholdN  = zeros(1,ncols);

  for iv = 1:ncols
    % if none survive the threshold is 0
    if maxrowID(iv); thresholdID(iv) = pValuesSorted(maxrowID(iv),iv); end
    if maxrowN(iv);  thresholdN(iv)  = pValuesSorted(maxrowN(iv),iv);  end
  end

  
function [] = testThis()

  % 100 voxels, progressively more active
  
  means{1} = zeros(1,100);
  means{2} = linspace(0,2,100);
  sigma{1} = eye(length(means{1}));
  sigma{2} = sigma{1};
  
  for c=1:2
    data{c} = mvnrnd(means{c},sigma{c},32);
  end

  pvalues = zeros(1,100);
  for v = 1:100
    [h,p] = ttest2(data{1}(:,v),data{2}(:,v),0.05,-1);
    pvalues(v) = p;
  end

  [signicant,threshold] = computeFDR( pvalues, 0.05 )
