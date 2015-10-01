function WholeBrain_MVPA(varargin)
  p = inputParser;
  p.KeepUnmatched = false;
  %% Parse and set parameters
  addParameter(p , 'debug'            , false     , @islogicallike );
  addParameter(p , 'SmallFootprint'   , false     , @islogicallike );
  addParameter(p , 'Gtype'            , []        , @ischar        );
  addParameter(p , 'normalize'        , false                      );
  addParameter(p , 'bias'             , false     , @islogicallike );
  addParameter(p , 'filters'          , []        , @ischarlike    );
  addParameter(p , 'filters_OR'       , []        , @ischarlike    );
  addParameter(p , 'target'           , []        , @ischar        );
  addParameter(p , 'data'             , []        , @ischarlike    );
  addParameter(p , 'data_var'         , 'X'       , @ischar        );
  addParameter(p , 'metadata'         , []        , @ischar        );
  addParameter(p , 'metadata_var'     , 'metadata', @ischar        );
  addParameter(p , 'cvfile'           , []        , @ischar        );
  addParameter(p , 'cv_var'           , 'CV'      , @ischar        );
  addParameter(p , 'cvscheme'         , []        , @isintegerlike );
  addParameter(p , 'cvholdout'        , []        , @isnumeric     );
  addParameter(p , 'orientation'      , []        , @ischar        );
  addParameter(p , 'diameter'         , []        , @isnumeric     );
  addParameter(p , 'overlap'          , []        , @isnumeric     );
  addParameter(p , 'shape'            , []        , @ischar        );
  addParameter(p , 'finalholdout'     , 0         , @isintegerlike );
  addParameter(p , 'lambda'           , []        , @isnumeric     );
  addParameter(p , 'alpha'            , []        , @isnumeric     );
  addParameter(p , 'AdlasOpts'        , struct()  , @isstruct      );
  addParameter(p , 'environment'      , 'condor'  , @ischar        );
  addParameter(p , 'SanityCheckData'  , []        , @ischar        );
  addParameter(p , 'COPY'             , []                         );
  addParameter(p , 'URLS'             , []                         );
  addParameter(p , 'executable'       , []                         );
  addParameter(p , 'wrapper'          , []                         );

  if nargin > 0
    parse(p, varargin{:});
  else
    try
      jdat = loadjson('./params.json');
    catch ME
      disp('Current directory does not contain "params.json".');
      rethrow(ME);
    end
    fields = fieldnames(jdat);
    jcell = [fields'; struct2cell(jdat)'];
    parse(p, jcell{:});
  end

  % private function.
  required = {'Gtype','data','metadata','cvfile','cvscheme','cvholdout','finalholdout'};
  assertRequiredParameters(p.Results,required);

  DEBUG            = p.Results.debug;
  SmallFootprint   = p.Results.SmallFootprint;
  Gtype            = p.Results.Gtype;
  normalize        = p.Results.normalize;
  BIAS             = p.Results.bias;
  filter_labels    = p.Results.filters;
  filter_labels_or = p.Results.filters_OR;
  target           = p.Results.target;
  datafile         = p.Results.data;
  data_var         = p.Results.data_var;
  cvfile           = p.Results.cvfile;
  cv_var           = p.Results.cv_var;
  cvscheme         = p.Results.cvscheme;
  cvholdout        = p.Results.cvholdout;
  orientation      = p.Results.orientation;
  diameter         = p.Results.diameter;
  overlap          = p.Results.overlap;
  shape            = p.Results.shape;
  finalholdoutInd  = p.Results.finalholdout;
  metafile         = p.Results.metadata;
  metadata_var     = p.Results.metadata_var;
  lambda           = p.Results.lambda;
  alpha            = p.Results.alpha;
  opts             = p.Results.AdlasOpts;
  environment      = p.Results.environment;
  SanityCheckData  = p.Results.SanityCheckData; %#ok<NASGU>

  % Check that the correct parameters are passed, given the desired algorithm
  [lambda, alpha] = verifyLambdaSetup(Gtype, lambda, alpha);

  % If values originated in a YAML file, and scientific notation is used, the
  % value may have been parsed as a string. Check and correct.
  if isfield(opts, 'tolInfeas')
    if ischar(opts.tolInfeas)
      opts.tolInfeas = sscanf(opts.tolInfeas, '%e');
    end
  end
  if isfield(opts, 'tolRelGap')
    if ischar(opts.tolRelGap)
      opts.tolRelGap = sscanf(opts.tolRelGap, '%e');
    end
  end

  % If cell array with one element, unpack element from cell.
  datafile = uncell(datafile);
  metafile = uncell(metafile);

  switch environment
  case 'condor'
    root = './';
    datadir = root;
    datafile = updateFilePath(datafile, datadir);
    metafile = updateFilePath(metafile, datadir);
    cvpath   = updateFilePath(cvfile,datadir);
    matfilename = 'results.mat';
    infofilename = 'info.mat';
  case 'chris'
    matfilename = 'results.mat';
    infofilename = 'info.mat';

  otherwise
    error('Environment %s not implemented.', environment);

  end

  %% Load metadata
  StagingContainer = load(metafile, metadata_var);
  metadata = StagingContainer.(metadata_var); clear StagingContainer;
  N = length(metadata);
  n = [metadata.nrow];
  d = [metadata.ncol];

  %% Compile filters
  rowfilter  = cell(N,1);
  colfilter  = cell(N,1);
  for i = 1:N
    if isempty(filter_labels)
      rowfilter{i} = true(1,n(i));
      colfilter{i} = true(1,d(i));
    else
      [rowfilter{i},colfilter{i}] = composeFilters(metadata(i).filter, filter_labels, 'logic', @all);
    end
  end

  rowfilter_or  = cell(N,1);
  colfilter_or  = cell(N,1);
  for i = 1:N
    if isempty(filter_labels_or)
      rowfilter_or{i} = true(1,n(i));
      colfilter_or{i} = true(1,d(i));
    else
      [rowfilter_or{i},colfilter_or{i}] = composeFilters(metadata(i).filter, filter_labels_or, 'logic', @any);
      rowfilter{i} = rowfilter{i} & rowfilter_or{i};
      colfilter{i} = colfilter{i} & colfilter_or{i};
    end
  end
  clear rowfilter_or colfilter_or

  %% Load CV indexes, and identify the final holdout set.
  % N.B. the final holdout set is excluded from the rowfilter.
  cvind = loadCV(cvfile, cv_var, cvscheme, rowfilter);
  for i = 1:N
    % Add the final holdout set to the rowfilter, so we don't even load
    % those data.
    finalholdout = cvind{i} == finalholdoutInd;
    rowfilter{i}(rowfilter{i}) = forceRowVec(rowfilter{i}(rowfilter{i})) & forceRowVec(~finalholdout);
    % Remove the final holdout set from the cvind, to match.
    cvind{i} = cvind{i}(~finalholdout);

    if finalholdoutInd > 0
      cvind{i}(cvind{i}>finalholdoutInd) = cvind{i}(cvind{i}>finalholdoutInd) - 1;
    end
  end
  if finalholdoutInd > 0
    % Adjust the cv holdout index(es) down if they are higher than the final holdout.
    if ~isempty(cvholdout)
      cvholdout(cvholdout>finalholdoutInd) = cvholdout(cvholdout>finalholdoutInd) - 1;
    end
  end

  %% Load data and select targets
  [X,subjix] = loadData(datafile, data_var, rowfilter, colfilter, metadata);
  metadata   = metadata(subjix);
  rowfilter  = rowfilter(subjix);
  colfilter  = colfilter(subjix);
  cvind      = cvind(subjix);

  Y = selectTargets(metadata, target, rowfilter);

  %% Include voxel for bias (conditional)
  fprintf('%-26s', 'Including Bias Unit');
  msg = 'NO';
  if BIAS
    msg = 'YES';
    X = addBiasUnit(X);
  end
  fprintf(': [%3s]\n', msg);

  % Normalize columns of X
  fprintf('%-26s', 'Normalizing columns of X');
  msg = 'NO';
  if normalize
    % This is handled later, after isolating the training set.
    msg = 'YES';
  end
  fprintf(': [%3s]\n', msg);

  % Final holdout index
  fprintf('%-26s', 'Final holdout index');
  fprintf(': [%3d]\n', finalholdoutInd);

  fprintf('Data loaded and processed.\n');

  %% Plug in the parameters and run
  switch Gtype
  case 'lasso'
    [results,info] = learn_category_encoding(Y, X, Gtype, ...
                      'lambda'         , lambda         , ...
                      'alpha'          , alpha          , ...
                      'cvind'          , cvind          , ...
                      'cvholdout'      , cvholdout      , ...
                      'normalize'      , normalize      , ...
                      'DEBUG'          , DEBUG          , ...
                      'SmallFootprint' , SmallFootprint , ...
                      'AdlasOpts'      , opts); %#ok<ASGLU>

  case 'searchlight'
    X = uncell(X);
    Y = uncell(Y)+1;
    cvind = uncell(cvind);
    colfilter = uncell(colfilter);

    % create a 3D binary mask
    z = strcmp({metadata.coords.orientation}, orientation);
    xyz = metadata.coords(z).xyz(colfilter,:);
    mask = coordsTo3dMask(xyz);

    % create the "meta" neighbourhood structure
    meta = createMetaFromMask(mask, 3);

    % translate X matrix into 3D+time
    classifier = 'gnb_searchmight';
    [am,pm] = computeInformationMap(X,Y,cvind,classifier,'searchlight', ...
                                meta.voxelsToNeighbours,meta.numberOfNeighbours,'testToUse','accuracyOneSided_permutation',1000);

    results.accuracy_map = am;
    results.pvalue_map = am;

  case 'soslasso'
    xyz = cell(numel(metadata),1);
    for ii = 1:numel(xyz)
      z = strcmp({metadata(ii).coords.orientation}, orientation);
      xyz{ii} = metadata(ii).coords(z).xyz(colfilter{ii},:);
    end
    G = coordGrouping(xyz, diameter, overlap, shape);
    [results,info] = learn_category_encoding(Y, X, Gtype, ...
                      'groups'         , G              , ...
                      'lambda'         , lambda         , ...
                      'alpha'          , alpha          , ...
                      'cvind'          , cvind          , ...
                      'cvholdout'      , cvholdout      , ...
                      'normalize'      , normalize      , ...
                      'DEBUG'          , DEBUG          , ...
                      'SmallFootprint' , SmallFootprint , ...
                      'AdlasOpts'      , opts); %#ok<ASGLU>
  end
  fprintf('Saving:\n');
  fprintf('\t%s\n',matfilename);
  fprintf('\t%s\n',infofilename);

  %% Revise cv indexes
  % Add the final holdout index to all results.
  [results.finalholdout] = deal(finalholdoutInd);
  % Adjust the cvholdout indexes to accomodate the final holdout index.
  if isfield(results,'cvholdout')
    cvholdout = [results.cvholdout];
    z = cvholdout >= finalholdoutInd;
    cvholdout(z) = cvholdout(z) + 1;
    cvholdout = mat2cell(cvholdout(:),ones(numel(cvholdout),1));
    [results.cvholdout] = deal(cvholdout{:});
  end

  %% Save results
  save(matfilename,'results');
  fprintf('Done!\n');
end

%% Local functions
function [lambda, alpha] = verifyLambdaSetup(Gtype, lambda, alpha)
% Each algorithm requires different lambda configurations. This private
% function ensures that everything has been properly specified.
  switch Gtype
  case 'lasso'
    if ~isempty(alpha)
      warning('Lasso does not use the alpha parameter. It is being ignored.');
    end
    assert(~isempty(lambda) , 'Lasso requires lambda.');
    alpha  = [];

  case 'iterativelasso'
    if ~isempty(alpha)
      warning('Lasso does not use the alpha parameter. It is being ignored.');
    end
    assert(~isempty(lambda) , 'Iterative Lasso requires lambda.');
    alpha  = [];

  case 'soslasso'
    assert(~isempty(lambda) , 'SOS Lasso requires lambda.');
    assert(~isempty(alpha)  , 'SOS Lasso requires alpha.');
  end
end

function assertRequiredParameters(params, required)
  N = length(required);
  for i = 1:N
    req = required{i};
    assert(isfield(params,req), '%s must exist in params structure! Exiting.');
    assert(~isempty(params.(req)), '%s must be set. Exiting.');
  end
end

function b = islogicallike(x)
  b = any(x == [1,0]);
end

function b = isintegerlike(x)
  b = mod(x,1) == 0;
end

function b = ischarlike(x)
  b = ischar(x) || iscellstr(x);
end

function r = rankind(ind)
  [~,ix] = sort(ind);
  [~,r]  = sort(ix);
end

function [X,subjix] = loadData(datafile,data_var,rowfilter,colfilter,metadata)
  % Load data for multiple subjects, and apply filters.
  datafile  = ascell(datafile);
  rowfilter = ascell(rowfilter);
  colfilter = ascell(colfilter);
  N         = length(datafile);
  subjix    = zeros(1,N);
  X         = cell(N,1);
  for i = 1:N
    fprintf('Loading %s from  %s...\n', data_var, datafile{i});
    subjid    = extractSubjectID(datafile{i});
    subjix(i) = find([metadata.subject] == subjid);
    tmp       = load(datafile{i}, data_var);
    X{i}      = tmp.(data_var); clear tmp;
    X{i}      = X{i}(rowfilter{subjix(i)},colfilter{subjix(i)});
  end
  X = uncell(X);
end

function Y = selectTargets(metadata, target, rowfilter)
  Y = {metadata.(target)};
  N = length(Y);
  for i = 1:N
    Y{i} = Y{i}(rowfilter{i});
  end
  Y = uncell(Y);
end

function cvind = loadCV(cvpath, cv_var, cvscheme, rowfilter)
  % In this dataset, trials are ordered the same way across subjects, so each
  % subject gets a copy of the same CV scheme.
  rowfilter = ascell(rowfilter);
  N = numel(rowfilter);
  tmp = load(cvpath, cv_var);
  CV = tmp.(cv_var); clear tmp;
  if N > 1
    tmp = {CV(:, cvscheme)};
    cvind = repmat(tmp, N, 1);
    for i=1:N
      cvind{i} = cvind{i}(rowfilter{i});
    end
  else
    rowfilter = uncell(rowfilter);
    cvind = CV(rowfilter, cvscheme);
  end
end

function id = extractSubjectID(datafile)
  [~,fname,~] = fileparts(datafile);
  str = regexp(fname,'[0-9]+','match');
  id = sscanf(str{1},'%d');
  if isempty(id)
    error('Failed to extract subject id from data filename %s. Exiting...', datafile);
  end
end

function filepath = updateFilePath(file, path, suffix)
  file     = ascell(file);
  N        = length(file);
  filepath = cell(1,N);
  for i = 1:N
    [~,f,e] = fileparts(file{i});
    if nargin>2
      f = sprintf('%s_%s%s',f,suffix,e);
    else
      f = sprintf('%s%s',f,e);
    end
    filepath{i} = fullfile(path,f);
  end
  if N == 1
    filepath = filepath{1};
  end
end
 
function X = shuffleData(X,rowfilter) %#ok<DEFNU>
  X = ascell(X);
  N = length(X);
  n = length(rowfilter{1});
  shuffledIndex = randperm(n);
  for i = 1:N
    ix   = rankind(shuffledIndex(rowfilter{i}));
    X{i} = X{i}(ix,:);
  end
  X = uncell(X);
end

function X = randomizeData(X) %#ok<DEFNU>
  X = ascell(X);
  N = length(X);
  for i = 1:N
    X{i} = randn(size(X{i}));
  end
  if N == 1
    X = X{i};
  end
  X = uncell(X);
end

function X = addBiasUnit(X)
  X = ascell(X);
  N = length(X);
  for i = 1:N
    X{i} = [X{i}, ones(size(X{i},1),1)];
  end
  X = uncell(X);
end

function C = ascell(X)
  if ~iscell(X)
    C = {X};
  else
    C = X;
  end
end

function M = uncell(C)
  % If C is a cell with one element, extract that element and descard the cell
  % casing.
  M = C;
  if iscell(C)
    if length(C) == 1
      M = C{1};
    end
  end
end

function b = isfunction_handle(x)
  if isa(x, 'function_handle');
    b = true;
  else
    b = false;
  end
end

function [rowfilter, colfilter] = composeFilters(filterset,labels,varargin)
  p = inputParser;
  p.KeepUnmatched = false;
  addParameter(p , 'logic', @all, @isfunction_handle);
  parse(p, varargin{:});
  combinelogic = p.Results.logic;

  labels = ascell(labels);

  % metadata.filter points to a structured array of filters.
  % First, force filters to a common orientation.
  for ii = 1:numel(filterset)
    filterset(ii).filter = forceRowVec(filterset(ii).filter);
  end

  % Then select the filters
  z = false(1,numel(filterset));
  for f = labels;
    z(strcmp(f, {filterset.label})) = true;
  end

  filters.row = filterset(z & [filterset.dimension]==1);
  filters.col = filterset(z & [filterset.dimension]==2);
  rowfilter = combinelogic(cat(1, filters.row.filter),1);
  colfilter = combinelogic(cat(1, filters.col.filter),1);
end

function r = forceRowVec(x)
  r = x(:)';
end

function mask = coordsTo3dMask(coords)
  dcoords = voxelspacing(coords, [3,3,3]); % [3,3,3] is just a seed.
  mcoords = min(coords);
  ijk = roundto(bsxfun(@minus, coords, mcoords),dcoords);
  ijk = roundto(bsxfun(@rdivide, ijk, dcoords),1)+1;
  ind = sub2ind(max(ijk), ijk(:,1), ijk(:,2), ijk(:,3));
  mask = false(max(ijk));
  mask(ind) = true;
end
