%% Introduction
%  ============
% To conceptually orient yourself, I highly recommend the first three
% chapters of this book by Trevor Hastie, Robert Tibshirani, and Martin
% Wainright
% (https://web.stanford.edu/~hastie/StatLearnSparsity_files/SLS_corrected_1.4.16.pdf).


% For example, imagine a study with 100 unique items, sampled equally from two
% categories or belonging to two experimental conditions. Further, imagine that
% there are 10,000 voxels in the cortex of this subject.
nitems = 100;
nvoxels = 10000;
X = randn(nitems, nvoxels);

% The assumption underlying sparse regression, including LASSO and all
% related methods like SOS Lasso, is the most of the voxels will be going
% unrelated things, and a small subset will be relevant. We can directly
% specify the contribution that each voxel will have to our target
% structure in a 10,000 element vector, b. Elements in b that are set to
% zero indicate that the voxel contributes nothing to the target structure.
b = zeros(10000, 1);
b(1:10) = 8;

% Given that we have defined our MRI data X and the contribution of each
% voxel to the structure, b, we can also generate a target structure, y.
y = X*b;