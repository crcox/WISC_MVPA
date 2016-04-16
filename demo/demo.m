% The data
% ========
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
nitems = 100;
nvoxels = 10000;
X = randn(nitems, nvoxels);

% Let y represent the category or condition labels of the items. This is the
% target structure that we will be modelling based on the data in X.
% Since all methods in this packaged are binary classifiers, y should be
% binary.
y= [true(50,1); false(50,1)];

% And, just for the sake of there being something to actually discover in our
% fake data, let's inject some structure into X. There will be 20 signal
% carrying voxels, and they will be clustered at the ``beginning'' of the
% dataset.
X(1:50,1:20) = X(1:50,1:20) + 2;

% Now, if we just wanted to fit a standard lasso, we could do something like:
[b,stats] = lassoglm(X, y,'binomial','CV',10);

% But, if we want to leverage WholeBrain_MVPA and do, say, SOS Lasso, there is
% some more setup we need to do.

% The metadata
% ============
% Targets
% -------
% Information about targets (i.e., possible y vectors) should be stored in a
% structure with 2 required fields: 'label' and 'targets'.
TARGETS(1) = struct('label','faces','targets',y);
TARGETS(2) = struct('label','places','targets',~y);

% Cross-validation
% ----------------
% The methods in WholeBrain_MVPA will not generate CV indexes for you. This is
% to help promote replicability of results. So you will need to provide a cross
% validation scheme ahead of time. If you are concerned about results being
% specific to a particular cross-validation scheme, you can specify multiple
% schemes.
%
% For example, let's set up 10 cross validation schemes.
nschemes = 10;
nfolds = 10;
SCHEMES = zeros(nitems, nschemes);
for iScheme = 1:nschemes
  c = cvpartition(y,'KFold', nfolds);
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
z = [true(200,1),false(9800,1)];
FILTERS(1) = struct('label','ROI01', 'dimension', 2, 'filter', z);
z = [true(98,1),false(2,1)];
FILTERS(2) = struct('label','GoodRows', 'dimension', 1, 'filter', z);

% Coordinates
% -----------
% It is useful (and in the case of SOS Lasso, essential) that the voxel
% coordinates be represented in the metadata structure. Like filters,
% coordinates are represented in a structure. The coordinate structure has 2
% required fields: 'orientation' and 'xyz'. The 'orientation' field is like the
% 'label' field in the filter structure and is used to look up particular
% coordinates. Labeling different coordinate spaces by 'orientation' is an AFNI
% convention. You don't have to use orientation codes like 'tlrc', 'orig', and
% 'mni', but that is my convention.
xyz = [1:nvoxels,1,1];
COORDS(1) = struct('orientation','orig','xyz',xyz);
COORDS(2) = struct('orientation','tlrc','xyz',xyz);

% Put it all together
% -------------------
% The metadata object compiles these three items, along with a couple other
% bits, into a single structure. The metadata structure has several required
% fields: 'subject', 'targets', 'filters', 'coords', 'cvind', 'nrow', 'ncol'.
% There will be a metadata structure for each subject, compiled into a
% structured array. Although in the example below subjects 100 and 101 are the
% same aside from their subject numbers, in practice they could be given
% different information.
metadata(1).subject = 100;
metadata(1).targets = TARGETS;
metadata(1).filters = FILTERS;
metadata(1).coords = COORDS;
metadata(1).cvind = SCHEMES;
metadata(1).nrow = nitems;
metadata(1).ncol = nvoxels;
metadata(2).subject = 101;
metadata(2).targets = TARGETS;
metadata(2).filters = FILTERS;
metadata(2).coords = COORDS;
metadata(2).cvind = SCHEMES;
metadata(2).nrow = nitems;
metadata(2).ncol = nvoxels;

% Save the data to disk
% =====================
% Despite having data and metadata organized properly in memory, before working
% with WholeBrain_MVPA we need to write the data to disk. The reason for this
% is that WholeBrain_MVPA is not written to be used interactively, but rather
% to facilitate to use in headless, batch applications particularly on
% distributed computing systems. WholeBrain_MVPA accepts paths to files on
% disk, as well as many other parameters.
% The data and metadata should be saved to a central location where it can be
% easily referenced.
% These files can be named whatever you like. You will be referencing them with
% explicit paths, and WholeBrain_MVPA does not make any assumptions about them.
% The program does assume that the *variable* names are X and metadata, but
% this default can be overwritten with certain parameters to WholeBrain_MVPA
% (data_var and metadata_var) if you prefer another convention.
subjects = [metadata.subject];
datadir = 'path/to/data'
for iSubj = 1:2
  s = subjects(iSubj);
  X = randn(nitems, nvoxels);
  X(1:50,1:20) = X(1:50,1:20) + 2;
  filename = sprintf('s%03d.mat', s);
  filepath = fullfile(datadir,filename);
  save(filepath, 'X');
end
save(fullfile(datadir,'metadata.mat'), 'metadata');

% Define a parameter file
% =======================
% WholeBrain_MVPA, despite being written as a Matlab function, is a very ugly
% function. First of all, it does not return anything. All results are written
% to disk. Likewise, although it is possible to invoke WholeBrain_MVPA from
% within a script or at the interactive terminal, it is designed to look for a
% parameter file if no arguments are provided. This is all makes
% WholeBrain_MVPA a bit counter-intuitive. However, these design choices make
% much more sense when considered in a distributed computing environment.
% WholeBrain_MVPA can be deployed to a system, along with a json file
% containing parameters, and it will parse the file and execute according to
% the instructions. It is designed to be executed with bare minimum interaction.
%
% Defining a parameter file is simple. See the documentation for a list of
% valid parameters. WholeBrain_MVPA reads json (http://www.json.org/), which is
% a widely used text-based syntax for representing structured data.
%
%           **The file must be named params.json**
%
% Below is an example parameter file, with valid json syntax. See the
% WholeBrain_MVPA docs for information about what each parameter does. This
% example will run SOS Lasso.
%
% {
%     "algorithm": "soslasso",
%     "bias": false,
%     "alpha": 0.4200556,
%     "lambda": 0.5863,
%     "shape": "sphere",
%     "diameter": 18,
%     "overlap": 9,
%     "cvscheme": 1,
%     "cvholdout": 0,
%     "finalholdout": 0,
%     "target": "faces"
%     "data": [
%         "/home/chris/data/FacePlaceObject/data/mat/handmade/TR5/jlp01.mat",
%         "/home/chris/data/FacePlaceObject/data/mat/handmade/TR5/jlp02.mat",
%         "/home/chris/data/FacePlaceObject/data/mat/handmade/TR5/jlp03.mat",
%         "/home/chris/data/FacePlaceObject/data/mat/handmade/TR5/jlp04.mat",
%         "/home/chris/data/FacePlaceObject/data/mat/handmade/TR5/jlp05.mat",
%         "/home/chris/data/FacePlaceObject/data/mat/handmade/TR5/jlp06.mat",
%         "/home/chris/data/FacePlaceObject/data/mat/handmade/TR5/jlp07.mat",
%         "/home/chris/data/FacePlaceObject/data/mat/handmade/TR5/jlp08.mat",
%         "/home/chris/data/FacePlaceObject/data/mat/handmade/TR5/jlp09.mat",
%         "/home/chris/data/FacePlaceObject/data/mat/handmade/TR5/jlp10.mat"
%     ],
%     "data_var": "X",
%     "normalize": "zscore",
%     "metadata": "/home/chris/data/FacePlaceObject/data/mat/handmade/TR5/metadata_TR5.mat",
%     "metadata_var": "metadata",
%     "orientation": "tlrc",
%     "filters": [
%         "rowfilter",
%         "colfilter"
%     ],
%     "SmallFootprint": false,
%     "debug": false,
% }
%
% Put the parameter file where you want to run the analysis. Paths can be
% relative with respect to where you execute WholeBrain_MVPA, but in most cases
% it will probably make sense for them to be absolute.

% Run WholeBrain_MVPA
% ===================
% With data and metadata structured properly and saved to disk, and with a
% parameter file named params.json in a folder where you would like to execute
% the analysis and return results, all that remains is to boot up Matlab in the
% directory that contains 'params.json' and execute WholeBrain_MVPA() at the
% command prompt. If you have compiled WholeBrain_MVPA into an executable (as
% would be necessary on a distributed computing cluster), you can execute
% Wholebrain_MVPA directly from the command line. In either case, it will read
% the parameter file and begin analysis. When it completes you will find a
% results.mat (or results.json) file in the directory where WholeBrain_MVPA was
% executed.

% Compile Results
% ===============
% If you are using WholeBrain_MVPA on a distributed computing cluster, you will
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
