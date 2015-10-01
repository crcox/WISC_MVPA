%
% compute similarity/distance between rows of two matrices
%
% in:
% - matrix A - #rowsA x #columns
% - matrix B - #rowsB x #columns (can be the same as A)
% - similarity measure - can be either 'correlation' (similarity) or 'euclidean' (distance)
%
%
% out:
% - similarity - #rowsA x #rowsB
%
% notes:
% - euclidean doesn't work yet
%
% history:
% - 2009 Nov 23 - created - fpereira@princeton.edu
%


function [S,indexMostSimilar,valueMostSimilar] = computeSimilarityByRow( varargin )

S = [];
if nargin < 3
  fprintf('syntax: computeSimilarityByRow( matrix A, matrix B (can be same as A), <''correlation''|''euclidean''>)\n');return;
end
matrixA = varargin{1};
matrixB = varargin{2};
measure = varargin{3};

[nrowsA,ncolsA] = size(matrixA);
[nrowsB,ncolsB] = size(matrixB);

if ncolsA ~= ncolsB
  fprintf('error: the two matrices must have the same number of columns\n');return;
else
  ncols = ncolsA;
end

switch measure
 case {'correlation'}
  [St,indexMostSimilar,valueMostSimilar] = rowiseSimilarity(matrixA',matrixB',0);
 case {'euclidean'}
  [St,indexMostSimilar,valueMostSimilar] = rowiseSimilarity(matrixA',matrixB',1);
 otherwise
  fprintf('error: the measure must be ''correlation'' or ''euclidean''\n');return;
end
S = St';
