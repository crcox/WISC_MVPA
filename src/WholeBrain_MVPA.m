function WholeBrain_MVPA(varargin)
  p = inputParser;
  p.KeepUnmatched = true;
  % ----------------------Set parameters-----------------------------------------------
  addParameter(p , 'debug'            , false     , @islogicallike );
  addParameter(p , 'SmallFootprint'   , false     , @islogicallike );
  addParameter(p , 'Gtype'            , []        , @ischar        );
  addParameter(p , 'normalize'        , false     , @islogicallike );
  addParameter(p , 'bias'             , false     , @islogicallike );
  addParameter(p , 'rowfilters'       , []        , @ischarlike    );
  addParameter(p , 'colfilters'       , []        , @ischarlike    );
  addParameter(p , 'target'           , []        , @ischar        );
  addParameter(p , 'data'             , []        , @ischarlike    );
  addParameter(p , 'dataVarname'      , 'X'       , @ischar        );
  addParameter(p , 'metadata'         , []        , @ischar        );
  addParameter(p , 'metaVarname'      , 'metadata', @ischar        );
  addParameter(p , 'cvfile'           , []        , @ischar        );
  addParameter(p , 'cvVarname'        , 'CV'      , @ischar        );
  addParameter(p , 'cvscheme'         , []        , @isintegerlike );
  addParameter(p , 'cvholdout'        , []        , @isintegerlike );
  addParameter(p , 'finalholdout'     , 0         , @isintegerlike );
  addParameter(p , 'lambda'           , []        , @isnumeric     );
  addParameter(p , 'alpha'            , []        , @isnumeric     );
  addParameter(p , 'AdlasOpts'        , struct()  , @isstruct      );
  addParameter(p , 'environment'      , 'condor'  , @ischar        );
  addParameter(p , 'SanityCheckData'  , []        , @ischar        );

  if nargin > 0
    parse(p, varargin{:});
  else
    jdat = loadjson('params.json');
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
  rowfilters       = p.Results.rowfilters;
  colfilters       = p.Results.colfilters;
  target           = p.Results.target;
  datafile         = p.Results.data;
  dataVarname      = p.Results.dataVarname;
  cvfile           = p.Results.cvfile;
  cvVarname        = p.Results.cvVarname;
  cvscheme         = p.Results.cvscheme;
  cvholdout        = p.Results.cvholdout;
  finalholdoutInd  = p.Results.finalholdout;
  metafile         = p.Results.metadata;
  metaVarname      = p.Results.metaVarname;
  lambda           = p.Results.lambda;
  alpha            = p.Results.lambda1;
  LambdaSeq        = p.Results.LambdaSeq;
  opts             = p.Results.AdlasOpts;
  environment      = p.Results.environment;
  SanityCheckData  = p.Results.SanityCheckData;

  % Check that the correct parameters are passed, given the desired algorithm
  [lam, alpha] = verifyLambdaSetup(Gtype, lambda, alpha);

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
    root = './';
    datadir = fullfile(root,'data');
    datafile = updateFilePath(datafile, datadir);
    metafile = updateFilePath(metafile, datadir);
    cvpath   = updateFilePath(cvfile,datadir);
    matfilename = 'results.mat';
    infofilename = 'info.mat';

  otherwise
    error('Environment %s not implemented.', environment);

  end

  % Load metadata
  tmp = load(metafile, metaVarname);
  metadata = tmp.(metaVarname); clear tmp;
  N = length(metadata);
  n = [metadata.nrow];
  d = [metadata.ncol];

  % Compile filters
  rowfilters = {metadata.rowfilters};
  colfilters = {metadata.rowfilters};
  rowfilter  = cell(N,1);
  colfilter  = cell(N,1);
  for i = 1:N
    rowfilter{i} = combineFilters(rowfilters, n);
    colfilter{i} = combineFilters(colfilters, d);
  end

  % Load CV indexes, and identify the final holdout set.
  % N.B. the final holdout set is excluded from the rowfilter.
  cvind = loadCV(cvfile, cvVarname, cvscheme, N);
  for i = 1:N
    finalholdout = cvind{i} == finalholdoutInd;
    rowfilter{i} = rowfilter{i} & ~finalholdout;
    cvind{i} = cvind{i}(rowfilter{i});
  end

  if finalholdoutInd > 0
    cvind(cvind>finalholdoutInd) = cvind(cvind>finalholdoutInd) - 1;
    % Adjust the cv holdout index(es) down if they are higher than the final holdout.
    if ~isempty(cvholdout)
      cvholdout(cvholdout>finalholdoutInd) = cvholdout(cvholdout>finalholdoutInd) - 1;
    end
  end

  % Load data and select targets
  [X,subjix] = loadData(datafile, dataVarname, rowfilter, colfilter);
  metadata   = metadata(subjix);
  rowfilter  = rowfilter(subjix);
  colfilter  = colfilter(subjix);

  Y = selectTargets(metadata, target, rowfilter);

  % metadata should contain:
  % - Filters (including outliers)
  % - Coordinates
  % NB: Use the input parameter ``sorted'' to indicate whether the data
  % provided are already sorted or not. Default is ``true''.
  % - Sort indexes (if the rows of data for each subject are in the same order,
  % these will just be 1:nrow. If each subject is in a different order for some
  % reason, then these indexes will reorder the rows into a common order.
  % - Unsort indexes (if a sort had to be applied, then these indexes should undo that sort).
  % - Warp data (for translating coordinates to Tlrc or MNI space).

  if ~isempty(SanityCheckData)
    disp('PSYCH! This is a simulation.');
    switch SanityCheckData
    case 'shuffle'
      disp('Shuffling rows of MRI data!!!')
      X = shuffleData(X, rowfilter);

    case 'random'
      disp('Generating totally random MRI data!!!')
      X = randomizeData(X);

    case 'use_shuffled'
      fprintf('Using pre-shuffled MRI data!!!\n')
      datafile = updateFilePath(datafile, [], 'shuffle');
      [X,subjix] = loadData(datafile, dataVarname, metadata);

    case 'use_random'
      fprintf('Using predefined random MRI data!!!\n', fname)
      datafile = updateFilePath(datafile, [], 'random');
      [X,subjix] = loadData(datafile, dataVarname, metadata);

    case 'real'
      disp('Using the true data, unaltered.');

    end
  end

  % Include voxel for bias
  fprintf('%-28s', 'Including Bias Unit:');
  msg = 'NO';
  if BIAS
    msg = 'YES';
    X = addBiasUnit(X);
  end
  fprintf('[%3s]\n', msg);

  % Normalize columns of X
  fprintf('%-28s', 'Normalizing columns of X:');
  msg = 'NO';
  if normalize
    msg = 'YES';
    % This is handled later, after isolating the training set.
    %%%%%
    %%%%%
    %%%%%

  end
  fprintf('[%3s]\n', msg);


  %% ----------------Visual, Audio or Semantic similarities and processing----------------
%  if ~isempty(SanityCheckTargets)
%    switch SanityCheckTargets
%    case 'shuffle'
%      disp('Shuffling target vector!!!')
%      simpath = fullfile(datadir,simfile);
%      allSimStructs = load(simpath);
%      S = allSimStructs.(simtype);
%      shidx = randperm(size(S,1));
%      S = S(shidx,shidx);
%
%    case 'use_shuffled'
%      fprintf('Using pre-shuffled target vector!!!')
%      [~,fname,~] = fileparts(simfile);
%      simfile = sprintf('%s_shuffled.mat',fname);
%      simpath = fullfile(datadir, simfile);
%      allSimStructs = load(simpath);
%      S = allSimStructs.(simtype);
%
%    case 'real'
%      disp('Using the true target vector, unaltered.')
%      simpath = fullfile(datadir,simfile);
%      allSimStructs = load(simpath);
%      S = allSimStructs.(simtype);
%    end
%  else
%    simpath = fullfile(datadir,simfile);
%    allSimStructs = load(simpath);
%    S = allSimStructs.(simtype);
%  end
%  S = S(rowfilter,rowfilter); clear allSimStructs;
%  S = S(~finalholdout, ~finalholdout);

  fprintf('Data loaded and processed.\n');

  %% ---------------------Setting algorithm parameters-------------------------
  [results,info] = learn_category_encoding(Y, X, Gtype, ...
                    'lambda'         , lambda         , ...
                    'alpha'          , alpha          , ...
                    'cvind'          , cvind          , ...
                    'cvholdout'      , cvholdout      , ...
                    'normalize'      , normalize      , ...
                    'DEBUG'          , DEBUG          , ...
                    'SmallFootprint' , SmallFootprint , ...
                    'AdlasOpts'      , opts); %#ok<ASGLU>

  fprintf('Saving:\n');
  fprintf('\t%s\n',matfilename);
  fprintf('\t%s\n',infofilename);

  save(matfilename,'-struct','results');
  save(infofilename,'-struct','info');

  fprintf('Done!\n');
end

function [lambda, alpha] = verifyLambdaSetup(Gtype, lambda, alpha);
% Each algorithm requires different lambda configurations. This private
% function ensures that everything has been properly specified.
  switch Gtype
  case 'lasso'
    if ~isempty(alpha)
      warning('Lasso does not use the alpha parameter. It is being ignored.');
    end
    assert(~isempty(lambda) , 'Lasso requires lambda.');
    lambda = lambda;
    alpha  = [];

  case 'iterativelasso'
    if ~isempty(alpha)
      warning('Lasso does not use the alpha parameter. It is being ignored.');
    end
    assert(~isempty(lambda) , 'Iterative Lasso requires lambda.');
    lambda = lambda;
    alpha  = [];

  case 'soslasso'
    assert(~isempty(lambda) , 'SOS Lasso requires lambda.');
    assert(~isempty(alpha)  , 'SOS Lasso requires alpha.');
    lambda = alpha;
    alpha  = alpha;
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

function filter = combineFilters(filters,n)
  if ~isempty(filters)
    filters = ascell(filters)
    N       = length(filters);
    filter  = false(n,1);
    for i = 1:N
      key = filters{i};
      z   = metadata.(key);
      filter(z) = true;
    end
  else
    filter = true(n,1);
  end
end

function r = rankind(ind)
  [~,ix] = sort(ind);
  [~,r]  = sort(ix);
end

function [X,subjix] = loadData(datafile,dataVarname,rowfilter,colfilter)
  datafile  = ascell(datafile);
  rowfilter = ascell(rowfilter);
  colfilter = ascell(colfilter);
  N         = length(datafile);
  subjix    = zeros(1,N);
  X         = cell(N,1);
  for i = 1:N
    fprintf('Loading %s from  %s...\n', dataVarname, datafile{i});
    subjid    = extractSubjectID(datafile{i});
    subjix(i) = find([metadata.subject] == subjid);
    tmp       = load(datafile{i}, dataVarname);
    X{i}      = tmp.(dataVarname); clear tmp;
    n         = size(X{i},1);
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

function cvind = loadCV(cvfile, cvVarname, cvscheme, N);
  tmp = load(cvpath, cvVarname);
  CV = tmp.(cvVarname); clear tmp;
  if N > 1
    tmp = {CV(:, cvscheme)};
    cvind = repmat(tmp, N, 1);
  else
    cvind = CV(:, cvscheme);
  end
end

function id = extractSubjectID(datafile)
  [path,fname,ext] = fileparts(datafile);
  id = sscanf(fname, 's%d');
  if ~isempty(id)
    error('Failed to extract subject id from data filename %s. Exiting...', datafile);
  end
end

function filepath = updateFilePath(file, path, suffix)
  file     = ascell(file);
  N        = length(file);
  filepath = cell(1,N);
  for i = 1:N
    [p,f,e] = fileparts(file{i});
    if nargin>2
      f = sprintf('%s_%s.%s',f,suffix,e);
    else
      f = sprintf('%s.%s',f,e);
    end
    filepath{i} = fullfile(path,);
  end
  X = uncell(X);
end

function X = shuffleData(X,rowfilter)
  X = ascell(X);
  N = length(X);
  n = length(rowfilter{1});
  shuffledIndex = randperm(n)
  for i = 1:N
    ix   = rankind(shuffledIndex(rowfilter{i}));
    X{i} = X{i}(ix,:)
  end
  else
  X = uncell(X);
end

function X = randomizeData(X)
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
  end
end

function M = uncell(C)
  % If C is a cell with one element, extract that element and descard the cell
  % casing.
  if iscell(C)
    if length(C) == 1
      M = C{1};
    end
  end
end
