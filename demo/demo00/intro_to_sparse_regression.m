%% Introduction to sparse regression
% Imagine a study that yields 100 patterns of neural activity,
% each composed of 1,000 voxels in the cortex of one subject. This can be
% described as a 100 x 1000 matrix, where each row is an example, and each
% column is a voxel.

X = randn(100, 1000);
save('demo_data_1000.mat', 'X');

load('demo_data_1000.mat', 'X');
[nitems, nvoxels] = size(X);

% The assumption underlying sparse regression, including Lasso and all
% related methods like SOS Lasso, is that most of the voxels will be
% contributing to unrelated processes, and a small subset will be relevant
% to what you are studying. We can directly specify the contribution that
% each voxel will have to our target structure in a 1,000 element vector,
% b. Elements in b that are set to zero indicate that the voxel contributes
% nothing to the target structure.

b = zeros(nvoxels, 1);
b(1:10) = 10;

% Given that we have defined our MRI data X and the contribution of each
% voxel to the structure, b, we can also generate a target structure, y.

y = X*b;

% b is our ground truth about the relationship between X and y. 

% Let's quickly split this into a training set of 90 examples and a holdout
% (test) set of 10 items. This will allow us to assess whether models are
% being overfit to the training set.

yt = y(1:90);
Xt = X(1:90,:);
yh = y(91:100);
Xh = X(91:100,:);

%% Standard regression, with perfect voxel pre-selection
% Using a standard linear model of the first 10 voxels, we can recover this
% relationship, training examples 1:90 and evaluating on examples 91:100.
%
% For the basic regression, we are going to use ordinary least squares
% regression (OLS) as opposed to generalized least squares (GLS). This is
% because GLS, in a sense, regularizes the data to help with
% multicoliniarity and other violated assumptions. It will only muddy the
% waters, for the sake of this demo, to use GLS. However, if you would like
% to repeat this demo with GLS regression, just change the following
% variable:

regression_type = 'OLS';

switch regression_type
    case 'OLS' % ordinary least squares
        b_lm = (Xt(:,1:10)'*Xt(:,1:10))\(Xt(:,1:10)'*yt);
        a0_lm = 0; % simplification: this model should not need bias.
    case 'GLS' % 
        b_lm = glmfit(Xt(:,1:10),yt);
        a0_lm = b_lm(1);
        b_lm(1) = [];
end
ytz = Xt(:,1:10) * b_lm(1:10) + a0_lm;
yhz = Xh(:,1:10) * b_lm(1:10) + a0_lm;
        
header = sprintf('Standard regression (%s), perfect pre-selection', regression_type);
DisplayModelSummary( header, b(1:10), b_lm(1:10), yt, ytz, yh, yhz );

% This shows the expected result: if we know exactly which voxels relate to
% the dependent variable, and we have a sufficient number of observations
% upon which the model can be fit, then the weights can be recovered. The
% whole point, though, is that in practice we do not know which of the
% 1,000 voxels relate to the dependent variable.

%% Standard regression, no voxel pre-selection
% If all 1,000 voxels are included in the analysis, then the model will be
% "overfit" to the training set. There is no longer a unique solution that
% fits the data, and the model can be fit to noise. This results in a small
% proportion of variance explained among the holdout items, and an
% estimated weight vector that does not at all resemble the true weight
% vector used to generate the target data.

switch regression_type
    case 'OLS' % ordinary least squares
        b_lm_overfit = (Xt'*Xt)\(Xt'*yt);
        a0_lm_overfit = 0; % simplification: this model should not need bias.
    case 'GLS' % 
        b_lm_overfit = glmfit(Xt,yt);
        a0_lm_overfit = b_lm(1);
        b_lm(1) = [];
end
ytz = Xt * b_lm_overfit + a0_lm_overfit;
yhz = Xh * b_lm_overfit + a0_lm_overfit;

header = sprintf('Standard regression (%s), all 1,000 voxels', regression_type);
DisplayModelSummary( header, b, b_lm_overfit, yt, ytz, yh, yhz );

% This model neither generalizes to the holdout set, nor identifies the
% appropriate voxels.

%% Ridge regression
% If voxel pre-selection is not desirable, then we can try to tackle the
% problem with regularized regression techniques. Regularized regression
% involves adding additional penalties to the the algorithm that finds the
% weights that define the relationship between the data and the model.
% Different penalties bias the weights in different ways. Each penalty
% attempts to overcome the problem of having lots of voxels and relatively
% few neural patterns to train on.
%
% Not all penalties will result in a sparse solution. A common example of a
% regularized regression technique that can improve generalization without
% providing a sparse solution is ridge regression.
%
% We'll use the ridge regression algorithm implemented as part of GLMNet,
% to help quickly pick a suitable hyperparameter. All regularized
% regression techniques will have one or more hyperparameters that scale
% the severity of the supplementary penalties that enforce the
% regularization.

% This will just add the GLMnet code to your path, if you have not already.

if ~exist('cvglmnet','file')
    addpath(GetFullPath(fullfile(pwd,'..','..','dependencies','glmnet')));
end

opts = glmnetSet();
opts_ridge_cv = opts;
opts_ridge_cv.alpha = 0;
fitobj_ridge_cv = cvglmnet(Xt,yt,'gaussian', opts_ridge_cv, 'mse', 10);

opts_ridge = glmnetSet(struct(...
    'alpha', 0, ...
    'lambda', fitobj_ridge_cv.lambda_min));

fitobj_ridge = glmnet(Xt,yt,'gaussian', opts_ridge);

a0_ridge = fitobj_ridge.a0; 
b_ridge = fitobj_ridge.beta;
ytz = glmnetPredict(fitobj_ridge, Xt);
yhz = glmnetPredict(fitobj_ridge, Xh);

header = sprintf('Ridge regression, all 1,000 voxels');
DisplayModelSummary( header, b, b_ridge, yt, ytz, yh, yhz );

% In this example, the model obtained via ridge regression generalizes
% somewhat better than the solution obtained by OLS regression, but still
% is only able to account for 12% of the variance. However, there is no way
% to tell which weights are contributing to this better generalization
% performance. It is not possible* to draw conclusions about which voxels
% are most important using ridge.

% * No simple off the shelf way, anyhow.

%% Lasso regression
% Sparse methods, while still biased estimates, will assign only a small
% number of non-zero weights when fitting a model. Since voxels that have
% zero weights are by definition not contributing to the model, we can at
% least say which voxels are not contributing to the model's ability to
% generalize (those with zero weights). Some of the selected voxels may be
% spuriously correlated noise that fit well to the training set and not to
% the holdout set, so we cannot say that all voxels that are selected are
% contributing to above chance classification.
%
% <https://arxiv.org/pdf/0809.2932.pdf Stability selection> and
% <https://arxiv.org/pdf/1403.4296.pdf permutation testing> are two methods
% which are being experimented with in the literature which may help
% separate the wheat from the chaff among the voxels selected by Lasso and
% other sparse methods.
%
% The more salient issue with Lasso and other sparse methods is that you
% cannot be sure that voxels that were not selected would not have also
% been predictive on the holdout set. Sparse methods, generally speaking,
% achieve sparsity by trying to do the most explanation with the fewest
% number of features (in our case, voxels). It follows, then, that if two
% or more features are highly correlated with each other, a sparse solution
% will select only one of them, since they contribute similar information.
% How to obtain sparse solutions that retain correlated features in the
% service of obtaining better estimates of the true weights is an active
% topic of research in optimization and machine learning. The ordered
% weighted Lasso (OWL) is one example, being pursued by researchers here at
% Wisconsin (including Rob Nowak and Urvashi Oswal, among others).
%
% We'll use the Lasso algorithm implemented as part of GLMNet,
% to help quickly pick a suitable hyperparameter. All regularized
% regression techniques will have one or more hyperparameters that scale
% the severity of the supplementary penalties that enforce the
% regularization.

% This will just add the GLMnet code to your path, if you have not already.

if ~exist('cvglmnet','file')
    addpath(GetFullPath(fullfile(pwd,'..','..','dependencies','glmnet')));
end

opts = glmnetSet();
opts_lasso_cv = opts;
opts_lasso_cv.alpha = 1;
fitobj_lasso_cv = cvglmnet(Xt,yt,'gaussian', opts_lasso_cv, 'mse', 10);

opts_lasso = glmnetSet(struct(...
    'alpha', 1, ...
    'lambda', fitobj_lasso_cv.lambda_min));

fitobj_lasso = glmnet(Xt,yt,'gaussian', opts_lasso);

a0_lasso = fitobj_lasso.a0; 
b_lasso = fitobj_lasso.beta;
ytz = glmnetPredict(fitobj_lasso, Xt);
yhz = glmnetPredict(fitobj_lasso, Xh);

header = sprintf('Lasso regression, all 1,000 voxels');
DisplayModelSummary( header, b, b_lasso, yt, ytz, yh, yhz );

%% Conclusion
% In this toy example with a high signal-to-noise ratio and where signal
% carrying voxels are uncorrelated to each other, Lasso performs
% essentially perfectly at recovering the true signal. This is not meant to
% inflate your confidence in Lasso in practice, and please use caution when
% interpretting sparse models fit to your true data!
%
% To conceptually orient yourself to regularized regression, I highly
% recommend the first three chapters of "Statistical Learning with
% Sparsity: The Lasso and Generalizations" by Trevor Hastie, Robert
% Tibshirani, and Martin Wainright
% (https://web.stanford.edu/~hastie/StatLearnSparsity_files/SLS_corrected_1.4.16.pdf).