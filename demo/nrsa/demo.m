% TUTORIAL: Run WholeBrain_RSA in a way suitable for distributed computing
% ========================================================================

%% The data
%  ========
% The fMRI data should be formatted in a time x voxel matrix, X. Each row of
% this matrix is a training example, and so should include all voxels that the
% model may be trained on/evaluated with. By ``time'', I also mean ``item''.
% That is, the method does not require that the data is the raw time series. In
% fact, in my experience it is more common to first peform an item-wise
% deconvolution which will result in a single volume of beta weights for each
% item. In that case, the item x voxel matrix will contain the fitted betas.

% For example, imagine a study with 100 unique items, sampled equally from two
% categories or belonging to two experimental conditions. Further, imagine that
% there are 10,000 voxels in the cortex of this subject.
nItems = 100;
nVoxels = 1000;
X = randn(nItems, nVoxels);

% Let S represent the target similarity structure. This might correspond to
% conditions or items, but the number of rows and columns in S must match the
% number of rows in X. This might involve averaging rows of X together.
Y = struct('visual',randn(100,3),'semantic',randn(100,8));
VISUAL = corr(Y.visual');
SEMANTIC = corr(Y.semantic');

%% The metadata
%  ============
% Targets
% -------
% Information about targets (i.e., possible y vectors) should be stored in a
% structure with 5 required fields: 'label', 'similarity', 'sim_source', 'sim_metric', and 'target'.
TARGETS(1).label = 'visual';
TARGETS(1).type = 'similarity'; % {'category','similarity'}
TARGETS(1).sim_source = 'chamfer';
TARGETS(1).sim_metric = 'chamfer';
TARGETS(1).target = VISUAL;

TARGETS(2).label = 'semantic';
TARGETS(2).type = 'similarity'; % {'category','similarity'}
TARGETS(2).sim_source = 'featurenorms';
TARGETS(2).sim_metric = 'cosine';
TARGETS(2).target = SEMANTIC;

% Cross-validation
% ----------------
% The methods in WholeBrain_RSA will not generate CV indexes for you. This is
% to help promote replicability of results. So you will need to provide a cross
% validation scheme ahead of time. If you are concerned about results being
% specific to a particular cross-validation scheme, you can specify multiple
% schemes.
%
% For example, let's set up 10 cross validation schemes.
nschemes = 10;
nfolds = 10;
SCHEMES = zeros(nItems, nschemes);
for iScheme = 1:nschemes
  c = cvpartition(nItems,'KFold', nfolds);
  for iFold = 1:nfolds
    SCHEMES(:,iScheme) = SCHEMES(:,iScheme) + (test(c, iFold) * iFold);
  end
end

% Filters
% -------
% You may want to be able to select/exclude subsets of voxels and items without
% needed to make multiple copies of the data. By specifying filters, you can
% pre-specify these subsets and apply them programmatically.
% A filter is represented as a structure with 3 required fields: label,
% dimension, and filter. The label names the filter so that it can be easily
% referenced later, dimension encodes whether the filter applies to rows (1) or
% columns (2) of X. The filter is a binary vectory that represents the filter
% itself.
% Here, lets set up (totally arbitrarily) a ROI filter and a filter to exclude
% outliers.
z = [true(500,1);false(9500,1)];
FILTERS(1) = struct('label','ROI01', 'dimension', 2, 'filter', z);
z = [true(98,1);false(2,1)];
FILTERS(2) = struct('label','GoodRows', 'dimension', 1, 'filter', z);

% Coordinates
% -----------
% It is useful (and in the case of searchlight, essential) that the voxel
% coordinates be represented in the metadata structure. Like filters,
% coordinates are represented in a structure. The coordinate structure has 2
% required fields: 'orientation' and 'xyz'. The 'orientation' field is like the
% 'label' field in the filter structure and is used to look up particular
% coordinates. Labeling different coordinate spaces by 'orientation' is an AFNI
% convention. You don't have to use orientation codes like 'tlrc', 'orig', and
% 'mni', but that is my convention.
xyz = [(1:nVoxels)',ones(nVoxels,1),ones(nVoxels,1)];
COORDS(1) = struct('orientation','orig','xyz',xyz);
COORDS(2) = struct('orientation','tlrc','xyz',xyz);

%% Put it all together
%  ===================
% The metadata object compiles these three items, along with a couple other
% bits, into a single structure. The metadata structure has several required
% fields: 'subject', 'targets', 'filters', 'coords', 'cvind', 'nrow', 'ncol',
% 'itemindex', and 'runindex'. The item index is used in case items are
% repeated and, for instance, may need to be averaged together. The run index
% is used to identify which block or scanner run each trial belongs to. There
% will be a metadata structure for each subject, compiled into a structured
% array. Although in the example below subjects 100 and 101 are the same aside
% from their subject numbers, in practice they could be given different
% information.
metadata(1).subject = 100;
metadata(1).targets = TARGETS;
metadata(1).filters = FILTERS;
metadata(1).coords = COORDS;
metadata(1).cvind = SCHEMES;
metadata(1).nrow = nItems;
metadata(1).ncol = nVoxels;
metadata(1).itemindex = 1:nItems;
metadata(1).runindex = ones(1,nItems);

metadata(2).subject = 101;
metadata(2).targets = TARGETS;
metadata(2).filters = FILTERS;
metadata(2).coords = COORDS;
metadata(2).cvind = SCHEMES;
metadata(2).nrow = nItems;
metadata(2).ncol = nVoxels;
metadata(2).itemindex = 1:nItems;
metadata(2).runindex = ones(1,nItems);

%% Save the data to disk
%  =====================
% Despite having data and metadata organized properly in memory, before working
% with WholeBrain_RSA we need to write the data to disk. The reason for this
% is that WholeBrain_RSA is not written to be used interactively, but rather
% to facilitate to use in headless, batch applications particularly on
% distributed computing systems. WholeBrain_RSA accepts paths to files on
% disk, as well as many other parameters.
% The data and metadata should be saved to a central location where it can be
% easily referenced.
% These files can be named whatever you like. You will be referencing them with
% explicit paths, and WholeBrain_RSA does not make any assumptions about them.
% The program does assume that the *variable* names are X and metadata, but
% this default can be overwritten with certain parameters to WholeBrain_RSA
% (data_var and metadata_var) if you prefer another convention.
subjects = [metadata.subject];
datadir = './shared';
if ~exist(datadir,'dir')
    mkdir(datadir);
end

% While writing data to file, we are going to also build some structure into
% the data that corresponds to the visual and semantic similarity structures
% generated above. Namely, I am going to take Y (the sqrt of S, resenting each
% item as a point in a d-dimensional space) and rotate the values around the
% origin (0,0,0).
Rotations = [30,60,90,120];
iRot = 0;
r = @(theta) [cos(theta), -sin(theta); sin(theta), cos(theta)];
for iSubj = 1:2
    s = subjects(iSubj);
    X = randn(nItems,nVoxels);
    % Rotate visual structure
    iRot = iRot + 1;
    theta = Rotations(iRot);
    Yv = Y.visual;
    n = size(Yv,2);
    R = eye(n);
    R(1:2, 1:2) = r(theta);
    Yv = Yv * R;
    % Rotate semantic structure
    iRot = iRot + 1;
    theta = Rotations(iRot);
    Ys = Y.semantic;
    n = size(Ys,2);
    R = eye(n);
    R(1:2, 1:2) = r(theta);
    Ys = Ys * R;
    % Add structure into X
    X(:,1:18) = X(:,1:18) + repmat(Yv,1,6);
    X(:,19:66) = X(:,19:66) + repmat(Ys,1,6);
    % Save data
    filename = sprintf('s%03d.mat', s);
    filepath = fullfile(datadir,filename);
    save(filepath, 'X');
end
save(fullfile(datadir,'metadata.mat'), 'metadata');

%% Define a parameter file
%  =======================
% WholeBrain_RSA, despite being written as a Matlab function, is a very ugly
% function. First of all, it does not return anything. All results are written
% to disk. Likewise, although it is possible to invoke WholeBrain_RSA from
% within a script or at the interactive terminal, it is designed to look for a
% parameter file if no arguments are provided. This is all makes
% WholeBrain_RSA a bit counter-intuitive. However, these design choices make
% much more sense when considered in a distributed computing environment.
% WholeBrain_RSA can be deployed to a system, along with a json file
% containing parameters, and it will parse the file and execute according to
% the instructions. It is designed to be executed with bare minimum interaction.
%
% Defining a parameter file is simple. See the documentation for a list of
% valid parameters. WholeBrain_RSA reads json (http://www.json.org/), which is
% a widely used text-based syntax for representing structured data.
%
%           **The file must be named params.json**
%
% To read and write json, you will need jsonlab
% (http://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files)
% which I have bundled with my code:
addpath('../../dependencies/jsonlab/');

% Put the parameter file where you want to run the analysis. Paths can be
% relative with respect to where you execute WholeBrain_RSA, but in most cases
% it will probably make sense for them to be absolute. The following should
% translate into a valid json file for the purpose of this demo.
params = struct('regularization','GrOWL2','bias',0,'tau',0.2,...
    'lambda1', 0.4200556, 'lambda', 0.5863, 'LambdaSeq', 'inf',...
    'cvscheme', 1,'cvholdout', 1:10, 'finalholdout', 0,...
    'target_label', 'visual', 'target_type','similarity', ...
    'sim_source', 'chamfer', 'sim_metric', 'chamfer',...
    'data', './shared/s100.mat', 'data_var', 'X',...
    'normalize_data', 'zscore', 'normalize_wrt', 'training_set', ...
    'normalize_target', 'center', 'metadata', './shared/metadata.mat',...
    'metadata_varname', 'metadata', 'orientation', 'tlrc',...
    'filters', {{'ROI01','GoodRows'}}, 'SmallFootprint', false,...
    'debug', false, 'SaveResultsAs','mat', 'subject_id_fmt', 's%03d.mat');

savejson('',params,'FileName','params.json','ForceRootName',false);

%% Run WholeBrain_RSA
%  ===================
% With data and metadata structured properly and saved to disk, and with a
% parameter file named params.json in a folder where you would like to execute
% the analysis and return results, all that remains is to boot up Matlab in the
% directory that contains 'params.json' and execute WholeBrain_RSA() at the
% command prompt. If you have compiled WholeBrain_RSA into an executable (as
% would be necessary on a distributed computing cluster), you can execute
% Wholebrain_MVPA directly from the command line. In either case, it will read
% the parameter file and begin analysis. When it completes you will find a
% results.mat (or results.json) file in the directory where WholeBrain_RSA was
% executed.
addpath('../../src/')
WholeBrain_RSA()


%% Run WholeBrain_RSA: Searchlight
%  ===============================
% Put the parameter file where you want to run the analysis. Paths can be
% relative with respect to where you execute WholeBrain_RSA, but in most cases
% it will probably make sense for them to be absolute. The following should
% translate into a valid json file for the purpose of this demo.

addpath('../../simitar-rsa')
params = struct('regularization', 'GrOWL2', 'bias', false,'tau',0.2,...
    'searchlight', true, 'slShape', 'cube', 'slSim_Measure', 'nrsa',...
    'slRadius', 8, 'slPermutationType', 'overMatrics', 'slPermutations', 1000,...
    'lambda1', 0.4200556, 'lambda', 0.5863, 'LambdaSeq', 'inf',...
    'cvscheme', 1,'cvholdout', 1:10, 'finalholdout', 0,...
    'target', 'visual', 'sim_source', 'chamfer', 'sim_metric', 'chamfer',...
    'data', './shared/s100.mat', 'data_var', 'X',...
    'normalize', 'stdev', 'metadata', './shared/metadata.mat',...
    'metadata_varname', 'metadata', 'orientation', 'tlrc',...
    'filters', {{'ROI01','GoodRows'}}, 'SmallFootprint', false,...
    'debug', false, 'SaveResultsAs','json');

savejson('',params,'FileName','params.json','ForceRootName',false);

WholeBrain_RSA()

%% Compile Results
%  ===============
% If you are using WholeBrain_RSA on a distributed computing cluster, you will
% quickly find that the volume of results is difficult to manage effectively. I
% have written some utility functions in Wholebrain_MVPA/util that attempt to
% facilitate common actions, like loading data from many jobs into a single
% matlab structure, writing tables of data, dumping coordinates of identified
% voxels, etc.
% Alternatively, you may find that your volume of data demands a database
% solution. Although the default is to return data in .mat files, which makes
% it easy to read back into matlab, results can also be output in json format
% which facilitates storing in a SQL or NoSQL database like MongoDB. Setting up
% such a database solution is far beyond the scope of this demo, but the squall
% project (github.com/ikinsella/squall) is a developing solution that utilizes
% MongoDB to great effect.
