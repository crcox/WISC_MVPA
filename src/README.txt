README file for the Sparse Overlapping Groups lasso test codes

This file contains code for performing Least Squares regression or classification under the Logistic loss using the Sparse Overlapping Group Lasso (SOGlasso), as described in the following papers:

1. Rao, Nikhil, et al. "Sparse overlapping sets lasso for multitask learning and its application to fmri analysis." Advances in Neural Information Processing Systems. 2013.

2. Rao, Nikhil, et al. "Classification with Sparse Overlapping Groups" arXiv preprint  (arXiv:1402.4512v2). 2014



_________________________________________________________________________________
USAGE
_________________________________________________________________________________

The main test code is test_soglasso.m . This file creates a toy signal with
overlapping groups, and obtains Gaussian measurements. Then, we use the
SOGlasso to recover the signal under both least squares and logistic regression
settings.
_________________________________________________________________________________
MAIN CODES:
_________________________________________________________________________________

SOG_least_squares.m

Solves the SOGlasso least squares problem. The inputs required are explained in
the help section (help SOG_least_squares in the MATLAB command prompt). At the
end, a final debasing step is performed to remove bias in the solution.


SOG_logistic.m

Identical to the least squares case, but for the logistic loss.

_________________________________________________________________________________
UTILITY CODES
_________________________________________________________________________________


replicate_matrix.m

performs the feature replication to convert the problem to a non overlapping
sparse group lasso. The other outputs are designed to make the main codes more
efficient.

multi_transpose.m

transposes a cell array (in the multi task learning setting).
