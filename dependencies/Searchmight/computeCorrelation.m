%
% compute correlation between columns of one or two matrices
%
% in:
% - matrix A - #points x #columnsA
%   optional:
%   - matrix B - #points x #columnsB (if absent, the second matrix is the same as the first)
%
% out:
% - correlation - #columnsA x #columnsB
%
% notes:
%
% history:
% - 2010 June 14 - created from prior code - fpereira@princeton.edu
%

function [S,indexMostSimilar,valueMostSimilar] = computeCorrelation( varargin )

S = [];
if nargin < 1
  fprintf('syntax: computeCorrelation( matrix A, [matrix B, defaults to being the same as A])\n');return;
end
matrixA = varargin{1};
sameMatrix = 0;
if nargin == 2
  matrixB = varargin{2};
else
  matrixB = matrixA;
  sameMatrix = 1;
end  

[nrowsA,ncolsA] = size(matrixA);
[nrowsB,ncolsB] = size(matrixB);

if nrowsA ~= nrowsB
  fprintf('error: the two matrices must have the same number of rows\n');return;
else
  nrows = nrowsA;
end

[St,indexMostSimilar,valueMostSimilar] = rowiseSimilarity(matrixA,matrixB,0,sameMatrix);
S = St'; clear St;



function [] = testThis()


% 1) compile the C code for use with your MATLAB

mex rowiseSimilarity.c

% 2) a test with 1000 "voxels" and 100 time points

X = randn(1000,1000);

tic; Squick = computeCorrelation(X); toc
tic; Scorr  = corrcoef(X); toc
tic; Scorr2 = squareform(pdist(X','correlation')); toc
Scorr2 = -1*(Scorr2 - 1);

Stmp = abs(Squick-Scorr2);
difference = sum(Stmp(:))


