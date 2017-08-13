%% Running SOS Lasso with *WholeBrain_MVPA*
% 

%% The data
% The fMRI data should be formatted so that each row is a training
% _example_, and the columns are the _features_ that express each example.
% An _example_ might correspond to the activation over voxels at a
% particular point in time, the beta values or t-values resulting from an
% initial univariate model of a design matrix convolved with an HRF, or
% anything else you like. Each feature, for our purposes, is a voxel. 
%
% * This matrix can be named whatever you like.
% * It must be saved in a .mat file, either on it's own or along with other
%   variables.
% * The .mat file itself can be named whatever you like, with the one
%   constraint that *WholeBrain_MVPA* will expect a subject number to be
%   present somewhere in the filename.
% * Numbers can be zero padded.
%
% In this demo, we'll call the matrix |X|, which will be the only variable we save to the .mat file.
%
% Imagine a study with 100 unique items, sampled equally from two
% categories or belonging to two experimental conditions. Further, imagine that
% there are 1,000 voxels in the cortex of this subject.
%
% Let |y| represent the category or condition labels of the items. This is
% the target structure that we will be modelling based on the data in |X|.
% Since all methods in this package are binary classifiers, |y| should be
% binary. We'll make |y| dependent on independent contributions from 10
% voxels.

nitems = 100;
nvoxels = 1000;
X = randn(nitems, nvoxels);

b = zeros(nvoxels, 1);
b(1:10) = 10;
y = X*b;
y = y > median(y);

tabulate(y);

%% Targets metadata.targets
% Information about targets (i.e., possible y vectors) should be stored in a
% structure with 5 required fields:
% 
% # |label|
% # |type|
% # |targets|
% # |sim_source|
% # |sim_metric|
%
% Only the first three are relevant for classification analyses, but all
% must be present.

TARGETS(1) = struct(...
  'label','faces',...
  'type','category',...
  'target',y, ...
  'sim_source',[],'sim_metric',[]);
TARGETS(2) = struct(...
  'label','places',...
  'type','category',...
  'target',~y, ...
  'sim_source',[],'sim_metric',[]);

%% Cross-validation
% Next, we'll define indexes for cross validation. A single
% cross-validation _scheme_ is a vector, containing the whole-numbers |1:k|,
% where |k| is the number of cross validation folds. 
% 
% The Matlab function |cvpartition|, which belongs to the Statistics Toolbox,
% is very helpful for generating cross-validation schemes. One of the major
% advantages of using |cvpartition| is that, if you give it a categorical
% target structure as its first argument, it will pick |k| holdout sets
% where each set has a balanced sample of each category.
%
% For example, try:
%
%   y = repmat((1:3)', 20, 1);
%   c = cvpartition(y, 'kfold', 10);
%   disp(c);
%   for i = 1:10
%       z = test(c, i);
%       disp('Holdout index:');
%       disp(find(z)');
%       tabulate(y(z));
%       fprintf('\n\n');
%   end
%
% You can pre-specify multiple cross-validation schemes. If you specify
% more than one, each individual scheme will be a column in a matrix. When
% defining your analysis, you will provide a |cvscheme|, which will simply
% be a column-index into this matrix you are defining.
%
% For example, let's set up 10 cross validation schemes, each defining a
% different data partition for 10-fold cross validation.
nschemes = 10;
nfolds = 10;
SCHEMES = zeros(nitems, nschemes);
for iScheme = 1:nschemes
    c = cvpartition(y,'KFold', nfolds);
    for iFold = 1:nfolds
        SCHEMES(:,iScheme) = SCHEMES(:,iScheme) + (test(c, iFold) * iFold);
    end
end

%% Filters
% You may want to be able to select/exclude subsets of voxels and items without
% needed to make multiple copies of the data. By specifying filters, you can
% pre-specify these subsets and apply them programmatically.
%
% A filter is represented as a structure with 3 required fields:
%
% # |label| - names the filter so that it can be easily referenced.
% # |dimension| - encodes whether the filter applies to rows (1) or
% columns (2) of X.
% # |filter| - a binary vectory that represents the filter
% itself.
%
% Here, lets set up (totally arbitrarily) two filters. The first will
% define a region of interest and the second will exclude outliers.
z = [true(500,1);false(9500,1)];
FILTERS(1) = struct('label','ROI01', 'dimension', 2, 'filter', z);
z = [true(98,1);false(2,1)];
FILTERS(2) = struct('label','GoodRows', 'dimension', 1, 'filter', z);

%% Coordinates
% Up to this point, we've considered our data strictly as a matrix, with
% rows as examples and columns as features. Of course, voxels exist in
% space, and in order to ultimately map weight vectors or information maps
% back into a brain-space for visualization, you'll need to know the
% coordinate of each voxel.
%
% Like filters, coordinates are represented in a structure that allows
% descriptive labeling, and also allows you to store the three
% index/coordinate types that AFNI will output (via *3dMaskDump*):
% 
% * |ind| - 1-dimensional index.
% * |ijk| - the integer-valued data-space coordinates.
% * |xyz| - the "real world" decimal valued coordinates (in _mm_).
% 
% Something a bit awkward and idiosyncratic about the coordinate structure
% is that the field that functions the same as |label| so many other places
% in the metadata structure is, here, called |orientation|. This is because
% I figured the only reason one would have multiple coordinate spaces
% per subject is to represent the result of multiple warps. Each coordinate
% space resulting from a warp, in AFNI parlance, is an orientation. So,
% that is the etymology of the |orientation| label.
%
% That all being said, I've been meaning to rename |orientation| to |label|
% for ages for the sake of internal consistency, and because
% *WholeBrain_MVPA* does not do anything special behind the scenes with the
% field. It is treated exactly as values in the |label| field are treated
% elsewhere. Which means there is nothing special about using _orig_ or
% _tlrc_ as values for the orientation field. (If there is part of the code
% the breaks if you don't use _orig_ or _tlrc_, that would qualify as a
% bug!)
% 

ind = (1:nvoxels)';
ijk = [(1:nvoxels)',ones(nvoxels,1),ones(nvoxels,1)];
xyz = bsxfun(@minus, ijk, [nvoxels/2, 1, 1]);

COORDS(1) = struct('orientation','orig','ind',ind,'ijk',ijk,'xyz',xyz);
COORDS(2) = struct('orientation','tlrc','ind',ind,'ijk',ijk,'xyz',xyz);

%% Putting it all together
% The metadata object compiles these three items, along with a couple other
% bits of information, into a single structure. The metadata structure has
% several required fields:
% 
% * |subject| - A numeric* subject ID.
% * |targets| - Which will contain something like the |TARGETS| structure
% defined above.
% * |filters| - Which will contain something like the |FILTERS| structure
% defined above.
% * |coords| - Which will contain something like the |COORDS| structure
% defined above.
% * |cvind| - Which will contain something like the |SCHEMES| matrix
% defined above.
% * |nrow| - The number of rows in the data matrix for this subject (before
% applying any of the filters contained in |metadata(s).filters|).
% * |ncol| - The number of columns in the data matrix for this subject (before
% applying any of the filters contained in |metadata(s).filters|).
% 
% There will be a metadata structure for each subject, compiled into a
% structured array. Although in the example below subjects 100 and 101 are
% the same aside from their subject numbers, in practice they could be
% given different information.

% Subject 100
metadata(1).subject = 100;
metadata(1).targets = TARGETS;
metadata(1).filters = FILTERS;
metadata(1).coords = COORDS;
metadata(1).cvind = SCHEMES;
metadata(1).nrow = nitems;
metadata(1).ncol = nvoxels;

% Subject 101
metadata(2).subject = 101;
metadata(2).targets = TARGETS;
metadata(2).filters = FILTERS;
metadata(2).coords = COORDS;
metadata(2).cvind = SCHEMES;
metadata(2).nrow = nitems;
metadata(2).ncol = nvoxels;

%% Save the data to disk
% Despite having data and metadata organized properly in memory, before working
% with *WholeBrain_MVPA* we need to write the data to disk. The reason for this
% is that *WholeBrain_MVPA* is not written to be used interactively, but rather
% to facilitate to use in headless, batch applications particularly on
% distributed computing systems. *WholeBrain_MVPA* accepts paths to files on
% disk, as well as many other parameters.
% The data and metadata should be saved to a central location where it can be
% easily referenced.
% These files can be named whatever you like. You will be referencing them with
% explicit paths, and *WholeBrain_MVPA* does not make any assumptions about them.
% The program does assume that the *variable* names are X and metadata, but
% this default can be overwritten with certain parameters to *WholeBrain_MVPA*
% (data_var and metadata_var) if you prefer another convention.
subjects = [metadata.subject];
datadir = './shared';
if ~exist(datadir,'dir')
    mkdir(datadir);
end
for iSubj = 1:2
  s = subjects(iSubj);
  X = randn(nitems, nvoxels);
  X(1:50,1:20) = X(1:50,1:20) + 2;
  filename = sprintf('s%03d.mat', s);
  filepath = fullfile(datadir,filename);
  save(filepath, 'X');
end
save(fullfile(datadir,'metadata.mat'), 'metadata');

%% Conclusion
% Setting up the metadata structure can be a bit of a hassle, and is by
% its nature labor intensive. It is important to take great care when
% setting it up, because the metadata structure is the primary
% data structure that *WholeBrain_MVPA* will reference when attempting to
% run analyses.
%
% The good news is, that you do not need to set it up anew for every
% analysis. On the contrary, in an ideal world you should only need to
% specify the metadata structure one time per project, unless new filters,
% targets, or coordinates become necessary, or the underlying data matrices
% themselves change in some way.
%
% Once the metadata structure is defined and saved to the hard-drive, we
% can get on with the more interesting work of specifying anayses.

%% Define a parameter file
% *WholeBrain_MVPA*, despite being written as a Matlab function, is a
% pretty atypical function.
%
% First of all, it does not return anything. All results are written
% to disk. In addition, while it is possible to invoke *WholeBrain_MVPA* from
% within a script or at the interactive terminal, it is designed to take instructions from a json-formatted parameter filelook for a
% parameter file if no arguments are provided. This all makes
% *WholeBrain_MVPA* a bit counter-intuitive.
% 
% However, these design choices make
% much more sense when considered in a distributed computing environment.
% *WholeBrain_MVPA* can be deployed to a system, along with a json file
% containing parameters, and it will parse the file and execute according to
% the instructions. It is designed to be executed with bare minimum interaction.
%
% Defining a parameter file is simple. See the documentation for a list of
% valid parameters. *WholeBrain_MVPA* reads json (<http://www.json.org/>), which is
% a widely used text-based syntax for representing structured data.
%
% *The file must be named params.json!*
%
% I call this out in bold because it is important... but in practice, it
% isn't something you will need to think much about. Another bit of code,
% part of my <https://github.com/crcox/condortools CondorTools> repository,
% called *setupJobs*, will write you params.json files for you. But we are
% not quite there yet.
%
% To read and write json, you will need jsonlab
% (<http://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files>)
% which is bundled with *WholeBrain_MVPA*:
if ~exist('savejson','file')
    addpath(GetFullPath(fullfile(pwd,'..','..','dependencies','jsonlab')));
end

% Put the parameter file where you want to run the analysis. Paths can be
% relative with respect to where you execute *WholeBrain_MVPA*, but in most cases
% it will probably make sense for them to be absolute. The following should
% translate into a valid json file for the purpose of this demo. 
params = struct('regularization', 'soslasso', 'bias', false, 'alpha', 0.1,...
    'lambda', 0.9, 'shape', 'sphere', 'diameter', 10, 'overlap', 5,...
    'cvscheme', 1,'cvholdout', 1:10, 'finalholdout', 0, 'target', 'faces',...
    'data', {{'./shared/s100.mat', './shared/s101.mat'}}, 'data_var', 'X',...
    'normalize', 'zscore', 'metadata', './shared/metadata.mat',...
    'metadata_var', 'metadata', 'orientation', 'tlrc', 'filters', ...
    {{'ROI01','GoodRows'}}, 'SmallFootprint', false, 'debug', false,...
    'SaveResultsAs','json','subject_id_fmt','s%d.mat');
savejson('',params,'FileName','params.json','ForceRootName',false);

%% Run *WholeBrain_MVPA*: SOS Lasso
%  ==============================
% With data and metadata structured properly and saved to disk, and with a
% parameter file named params.json in a folder where you would like to execute
% the analysis and return results, all that remains is to boot up Matlab in the
% directory that contains 'params.json' and execute _WholeBrain_MVPA()_ at the
% command prompt. If you have compiled *WholeBrain_MVPA* into an executable (as
% would be necessary on a distributed computing cluster), you can execute
% *WholeBrain_MVPA* directly from the command line. In either case, it will read
% the parameter file and begin analysis. When it completes you will find a
% results.mat (or results.json) file in the directory where *WholeBrain_MVPA* was
% executed.
if ~exist('WholeBrain_MVPA','file')
    addpath(GetFullPath(fullfile(pwd,'..','..','src')));
end
WholeBrain_MVPA()

%%  Run *WholeBrain_MVPA*: Searchlight
%  ================================
% Put the parameter file where you want to run the analysis. Paths can be
% relative with respect to where you execute WholeBrain_MVPA, but in most cases
% it will probably make sense for them to be absolute. The following should
% translate into a valid json file for the purpose of this demo. 
% params = struct('algorithm', 'soslasso', 'bias', false, 'alpha', 0.4200556,...
%     'lambda', 0.5863, 'shape', 'sphere', 'diameter', 18, 'overlap', 9,...
%     'cvscheme', 1,'cvholdout', 1:10, 'finalholdout', 0, 'target', 'faces',...
%     'data', {{'./shared/s100.mat', './shared/s101.mat'}}, 'data_var', 'X',...
%     'normalize', 'zscore', 'metadata', './shared/metadata.mat',...
%     'metadata_var', 'metadata', 'orientation', 'tlrc', 'filters', ...
%     {{'ROI01','GoodRows'}}, 'SmallFootprint', false, 'debug', false,...
%     'SaveResultsAs','json','subject_id_fmt','s%d.mat');
% savejson('',params,'FileName','params.json','ForceRootName',false);

% Compile Results
% ===============
% If you are using *WholeBrain_MVPA* on a distributed computing cluster, you will
% quickly find that the volume of results is difficult to manage effectively. I
% have written some utility functions in *WholeBrain_MVPA*/util that attempt to
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
