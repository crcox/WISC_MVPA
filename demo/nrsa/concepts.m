% OVERVIEW OF KEY CONCEPTS AND MOTIVATION
% =======================================
% First, let's take a very simple example with a 25 voxel brain. There are
% several important concepts at play. Let's imagine a study where participants
% made semantic judgements about 50 items. We have an estimate of the semantic
% similarity structure derived from, e.g., a set of feature norms. Each item is
% associated with, say, 300 features that can be either present or absent. The
% similarity among these feature norms is summarized by their pairwise euclidean
% distances, which will yield a 50x50 symmetric matrix with zeros on the diagonal.

% That is the scenario you should have in mind, but for the purposes of the
% simulated example we will actually generate the data so there is a known
% relationship between the ``MRI data'' (X) and the similarity structure (S) by
% generating the S from X via a set of weights, U.

% Rather than being a vector with lenght nVoxels, U is a nVoxel x nDimension
% matrix. Dimension here is related to the ``rank'' of S, where structured
% solutions will have lower rank. Put another way, in a singular value
% decomposition (e.g. principle components analysis) of a matrix, the matrix
% will be decomposed into as many components as their are features, but only a
% small number of them will have large eigenvalues. When you plot the
% eigen-values of a non-random matrix, you will typically see that the
% magnitude of the eigen-values drops precipitously after the first few. This
% means that most of the variance in the matrix can be explained in a few
% dimensions, which is another way of saying there is ``similarity structure''
% expressed within the matrix.

nItems = 50;
nVoxels = 25;
nDimensions = 3;
nSignificantVoxels = 25;

% U will encode the relationship between X and sqrt(S)=Y. This is a little
% trick. Y=sqrt(S), and U=sqrt(W). The relationship between S and X would be
% modeled as S=XWX', but this is a difficult problem to solve. Y=XU is much
% easier to solve, and is related via a sqrt to the solution we want and so
% this is how we do it in practice.
X = randn(nItems,nVoxels);
U = zeros(nVoxels, nDimensions);
U(1:nSignificantVoxels,:) = randn(nSignificantVoxels, nDimensions);
W = U*U';

% N.B. Y will have as many columns as U, and the number of columns corresponds directly to the ``rank'' of S, because S=YY'
Y = X * U;
S = Y*Y';

% Now we have a ground truth for S, and W, and we have X. In practice, we will
% have only S and X and we want to solve for W. Because we have more items than
% voxels in this toy example, we can solve this relatively simply.

% The following closed form solution for the weight estimates, Uz, is based on
% the matrix formalization of regression. Note that the `z' suffix will mark
% estimated values. The `1' suffix means that the predictions or values are
% about the un-trained holdout set. We will fit these parameters, holding out
% 20% of the items for a test set.
c = cvpartition(50, 'HoldOut', .2)
X1 = X(test(c),:);
Y1 = Y(test(c),:);
X2 = X(training(c),:);
Y2 = Y(training(c),:);
Uz = pinv(X1) * Y1;

% The following is one possible way to compute an error metric. It's important
% to note that this error metric, in this context, does not have a natural
% value associated with ``chance'' or zero-structure, and it should be
% determined (in practice) by permutation test. Note also that, in practice, we
% should have done cross validation.
S1 = Y1 * Y1';
Yz1 = X(test(c),:) * Uz;
Sz1 = Yz1 * Yz1';
err_Y = norm(Y1 - Yz1,'fro')/norm(Y1,'fro');
err = norm(S1 - Sz1,'fro')/norm(S1,'fro');

% err_Y and err are fundamentally related, but both are computed for
% completeness. Error is essentially zero.

% You might be wondering how this differs from a standard RSA, where the
% similarity structure in the MRI data is no modeled, but assessed directly. We
% can make this comparison explicit by simply computing the euclidean distance
% among the rows of X, and seeing who closely this structure corresponds to the
% true S. Here, I am computing a ``similarityStructureScore'', which is simply
% the sum over the element-wise products between the two matrices. Larger
% scores are better---as the similarity between the two matrices increases, so
% does the score. Imporant to note, again, that this is on a totally arbitrary
% scale, and so can only be interpretted relative to something else, or to
% permutation tests.
d1 = pdist(X1);
D1 = squareform(d1);
similarityStructureScore = zeros(1,2);
similarityStructureScore(1) = S1(:)' * D1(:);
similarityStructureScore(2) = S1(:)' * Sz1(:);

% It is clear that the similarity score between S and Sz is much, *much* higher
% than S and D. What is going on? In a standard RSA, there is no model. It
% assumes that the similarity structure is expressed equally over all voxels
% considered, and that all voxels are equally relevant to all underlying
% dimensions of the similarity space. Network RSA allows the relevance of each
% voxel to each dimension to be weighted, and the predicted similarity
% structure is related to the weighted sum of neural activity. By anology,
% consider multivoxel (pattern) classification analysis. In MVPA, the voxel
% activity is modeled to make category predictions, and does not assume that
% all voxels are equally important. Network RSA brings this flexibility to RSA.

% Another way of driving this point home is to realize that the standard RSA
% analysis is similar to assuming that the weight on all voxels is 1. With this
% in mind, we can "make predictions" using the "standard model", which is
% simply the identity matrix.
Si1 = X1*eye(nVoxels)*X1';
err = zeros(1,2);
err(1) = norm(S1 - Si1,'fro')/norm(S1,'fro');
err(2) = norm(S1 - Sz1,'fro')/norm(S1,'fro');

% This reiterates the comparison above: Standard RSA suggests there is no
% structure in the data, but Network RSA shows that the data actually express
% the structure nearly perfectly.

% At this point, many of the core concepts have been laid out and compared with
% more basic RSA approaches. Next, we consider another key aspect of Network RSA,
% namely that sparse signal can be identified among a large set of features.
% Put another way, this means that Network RSA can discover a subset of voxels
% within a large volume (perhaps all of cortex) that encode similarity
% structure.

nItems = 50;
nVoxels = 1000;
nDimensions = 3;
nSignificantVoxels = 25;

X = randn(nItems,nVoxels);
U = zeros(nVoxels, nDimensions);
U(1:nSignificantVoxels,:) = randn(nSignificantVoxels, nDimensions) * 2;
W = U * U';
Y = X * U;
S = Y * Y';

c = cvpartition(50, 'HoldOut', .2)
X1 = X(test(c),:);
Y1 = Y(test(c),:);
S1 = Y1 * Y1';
X2 = X(training(c),:);
Y2 = Y(training(c),:);
Uz = pinv(X2) * Y2;

Yz1 = X1 * Uz;
Sz1 = Yz1 * Yz1';
err_Y = norm(Y1 - Yz1,'fro')/norm(Y1,'fro');
err = norm(S1 - Sz1,'fro')/norm(S1,'fro');

% Now, the error is very high. How high? We would need to do a permutation test
% to have a point of reference.
nPerm = 1000;
perm_err = zeros(nPerm,1);
for iPerm = 1:nPerm
  Uz = pinv(X2(randperm(c.TrainSize),:)) * Y2;
  Yz1 = X1 * Uz;
  Sz1 = Yz1 * Yz1';
  perm_err(iPerm) = norm(S1 - Sz1,'fro')/norm(S1,'fro');
end
fprintf('Proportion of permutations with error less than CV error: %.3f\n',nnz(perm_err < err)/nPerm)

% Clearly, we are unable to do better than chance when the signal is contained
% in a sparse subset of voxels. The problem is that there is no unique
% solution, and the training set is being over-fit. One solution to the
% overfitting problem is LASSO. In this case, because there are multiple
% columns of W, the solution is actually to use Group LASSO. Group Lasso is
% implemented within the WholeBrain_RSA toolbox.
addpath('../src');

% LASSO involves a free parameter that scales how much the model is penalized
% for having non-zero weights. If you want matlab to find a good value for
% this, run the line before. Otherwise, you can just try a few values and see
% if one works. In production you'll want to do some sort of principled
% parameter search, but for now something in the 0--10 range will probably do.
[lam1, err_L1L2] = fminbnd(@(x) optimizeAdlas1(X,Y,c,x), 0, 10)

[Uz, info] = Adlas1(X2, Y2, lam1);

% This error looks better, but we'll need to do a new permutation test to be sure.
% This is a slower procedure, so I'm going to dial back the permutation count...
nPerm = 10;
perm_err_L1L2 = zeros(nPerm,1);
for iPerm = 1:nPerm
  disp(iPerm)
  [Uz, info] = Adlas1(X2, Y2, lam1);
  Uz = pinv(X2(randperm(c.TrainSize),:)) * Y2;
  Yz1 = X1 * Uz;
  Sz1 = Yz1 * Yz1';
  perm_err_L1L2(iPerm) = norm(S1 - Sz1,'fro')/norm(S1,'fro');
end
fprintf('Proportion of permutations with error less than CV error: %.3f\n',nnz(perm_err_L1L2 < err_L1L2)/nPerm)

% This permutation test suggests that the difference between the predicted and
% true similarity structure is smaller than would be expected by chance.

% One limitation of LASSO (and Group LASSO) is that is seeks the sparsest
% possible solution, and so if there are multiple correlated features, LASSO
% will select one and drop the others. Ordered Weighted LASSO (OWL) is an
% alternative that helps retain correlated voxels and so makes more sense in a
% neuroscience context. Group OWL (GrOWL) is implemented in the WholeBrain_RSA
% package. Of course, these simulated data are generated at random and so we
% shouldn't expect an increase, and perhaps may even see a decrease, in
% preformance on these i.i.d. simulated data.

% GrOWL (implemented in Adlas2) takes 2 lambda arguments. The first affects how
% much to care about correlated voxels, and the second effects how much to care
% about sparsity.
[lamvals, err_growl] = fminsearch(@(x) optimizeAdlas2(X,Y,c,x), [.001,1])

[Uz, info] = Adlas2(X2, Y2, lamvals(1), lamvals(2));

% Probably due to the i.i.d. nature of the data, this yields a high error.
% Again, though, we will not have context for this error term without a
% permutation test.
nPerm = 10;
perm_err_growl = zeros(nPerm,1);
for iPerm = 1:nPerm
  disp(iPerm)
  [Uz, info] = Adlas2(X2, Y2, lam1);
  Uz = pinv(X2(randperm(c.TrainSize),:)) * Y2;
  Yz1 = X1 * Uz;
  Sz1 = Yz1 * Yz1';
  perm_err_growl(iPerm) = norm(S1 - Sz1,'fro')/norm(S1,'fro');
end
fprintf('Proportion of permutations with error less than CV error: %.3f\n',nnz(perm_err_L1L2 < err_L1L2)/nPerm)

% We have now introduced Network RSA, and how it can be achieved with Group
% LASSO and GrOWL using functions in the WholeBrain_RSA package. However, the
% computational intensity of working with these methods may have impressed
% themselves upon you at this point. Parameter searching, k-fold cross
% validation, and permutation testing, mean looping over relatively long
% running processes thousands of times. For this reason, WholeBrain_RSA was
% written with distribued computing in mind.
