%% Setup coordinates
% Here we use a helper function to define coodinates for a realistically
% sized brain space, and then select a sphere of them centered on the
% center of the space.
xyz = ndcoord(1:30,1:30,1:60);
xyz = bsxfun(@minus, xyz, mean(xyz));
d = pdist2(xyz,[0,0,0]);
xyz = xyz(d<7,:); % adjust here to grow the "brain".
clear d;

%% Divide into subjects
% Use the brain coordinates designated above to generate some subjects.
nsubj = 5;
XYZ = cell(1, nsubj);
for j = 1:nsubj
  idx = randperm(size(xyz,1),floor(size(xyz,1)/3));
  XYZ{j} = xyz(idx,:);
end

%% Generate groups
% This function will panel the space with either cubes or spheres, with
% given diameter and amount of overlap. The diameter and overlap are
% defined in terms of the coordinate units (e.g., millimeters in the fMRI
% context). If overlap = 4, then the radii of two adjacent groups will
% overlap by at most 4 millimeters.
% Note that G is simply a matrix of cells, one row per group and one column
% per subject. Within each cell are indexes that point to columns in X,
% asigning those columns to that group for that subject. These groups could
% be constructed in any which way, but must be presented to SOS Lasso in
% this way. coordGrouping makes it easy to tile the space with groups, but
% they could be constructed manually, or based on other data.
diameter = 9;
overlap = 4.5;
shape = 'sphere';
G = coordGrouping(XYZ,diameter,overlap,shape);

%% Generate some data
X = cell(size(XYZ));
Y = cell(size(XYZ));
B = cell(size(XYZ));
nsig = 200;
p = 500;
ag = 1:10;
snr = 2;
for j = 1:nsubj
  d = size(XYZ{j},1);
  x = randn(p, d);
  aix = unique(cell2mat(G(ag,j)));
  aix = aix(randperm(length(aix),nsig));
  x(:,aix) = x(:,aix) * snr;
  b = randn(nsig,1) * snr;
  y = x(:,aix) * b;
  bb = zeros(d,1);
  bb(aix) = b;
  X{j} = x;
  Y{j} = sign(y);
  B{j} = bb;
  clear x y b bb
end

%% Run SOS Lasso
alpha = 0.6;
lambda  = 0.4;

tic;
[W, obj] = SOS_logistic(X, Y, alpha, lambda, G);
toc

%% Debias the solution
% The weights obtained through regularized regression are biased estimates
% of the true parameters. One way to counteract this is to use the
% regularized regression purely as feature selection, and then fit a
% an unbiased model using only those selected features.
Wopt = cell(size(W));
Yz = zeros(p,size(W,2));
for j = 1:size(W,2)
  d = length(W{j});
  Wopt{j} = zeros(d,1);
  z = W{j}~=0;
  x = X{j};
  x = x(:, z);
  Wopt{j}(z) = x\Y{j};
  Yz(:,j) = X{j} * Wopt{j};
end
