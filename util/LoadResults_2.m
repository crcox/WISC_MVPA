function [results, params] = LoadResults(varargin)
  p = inputParser();
  addParameter(p,'ResultDir','.',@ischar);
  addParameter(p,'DataDir','~/MCW/WholeBrain_RSA/data/avg',@ischar)
  addParameter(p,'orientation','mni')
  addParameter(p,'CVSchemesFile','CV_schemes.mat',@ischar)
  addParameter(p,'cv_var','CV',@ischar)
  addParameter(p,'MetadataFile','metadata.mat',@ischar)
  addParameter(p,'meta_var','metadata',@ischar)
  addParameter(p,'Filters',{'TrueFaces'},@iscell)
  addParameter(p,'ResultFile','',@ischar);
  addParameter(p,'ParamFile','',@ischar);
  addParameter(p,'SortJobs',false,@islogical)
  addParameter(p,'SkipFields',[])
  addParameter(p,'Query',[])
  parse(p,varargin{:});

  resultdir     = p.Results.ResultDir;
  datadir       = p.Results.DataDir;
  orientation   = p.Results.orientation;;
  cvsfile       = p.Results.CVSchemesFile;
  cv_var        = p.Results.cv_var;
  metafile      = p.Results.MetadataFile;
  meta_varname  = p.Results.meta_var;
  filter_labels = p.Results.Filters;
  sortjobs      = p.Results.SortJobs;
  resultfile    = p.Results.ResultFile;
  paramfile     = p.Results.ParamFile;
  SKIP          = p.Results.SkipFields;
  Queries       = p.Results.Query;
  allfiles      = dir(resultdir);
  alldirs       = allfiles([allfiles.isdir]);
  jobdirs       = SelectJobDirs(alldirs, 'root',resultdir, 'sort', sortjobs);

  if ~isempty(SKIP) % not used
    SkipStr = SKIP;
    SkipStr{end} = ['and ', SKIP{end}];
    SkipStr = strjoin(SkipStr, ', ');
    fprintf('Fields %s will be skipped.\n', SkipStr);
  end

  %% Load Metadata
  metapath = fullfile(datadir, metafile);
  StagingContainer = load(metapath, meta_varname);
  metadata = StagingContainer.(meta_varname);
  N = numel(metadata);

  %% Preallocate result structure. Assumption is that all result files are
  %identical size with same fields.
  i = 0;
  FileExists = false;
  while ~FileExists
    i = i + 1;
    jobdir     = fullfile(resultdir, jobdirs(i).name);
    resultpath = fullfile(jobdir, resultfile);
    FileExists = exist(resultpath, 'file') == 2;
  end
  tmp = load(resultpath);
  R = tmp.results;
  N = numel(jobdirs) * numel(R);
  results = R(1);
  results.job = 0;
  results(N).job = 0;
  disp(results)

  parampath   = fullfile(jobdir,paramfile);
  P = loadjson(parampath);
  Overwrite = struct('data',0,'cvholdout',0,'lambda',0,'alpha',0);
  if ~isempty(Queries) && ~isfield(Queries, 'data')
    Queries.data = P.data;
    Overwrite.data = 1;
  end
  if ~isempty(Queries) && ~isfield(Queries, 'cvholdout')
    Queries.cvholdout = P.cvholdout;
    Overwrite.cvholdout = 1;
  end
  if ~isempty(Queries) && ~isfield(Queries, 'lambda')
    Queries.lambda = P.lambda;
    Overwrite.lambda = 1;
  end
  if ~isempty(Queries) && ~isfield(Queries, 'alpha')
    Queries.alpha = P.alpha;
    Overwrite.alpha = 1;
  end

  %% Loop over job dirs
  n = length(jobdirs);
  nchar = 0;
  fprintf('Loading job ');
  iii = 0;
  iiir = 0;
  for i = 1:n;
    fprintf(repmat('\b', 1, nchar));
    nchar = fprintf('%d of %d', i, n);

    %% load parameter file
    jobdir      = fullfile(resultdir, jobdirs(i).name);
    if nargout > 1
      parampath   = fullfile(jobdir,paramfile);
      tmp         = loadjson(parampath);
      tmp.jobdir  = jobdir;
      params(i)   = tmp;
      clear tmp;
      P = params(i);
      if Overwrite.data
        Queries.data = P.data;
      end
      if Overwrite.cvholdout
        Queries.cvholdout = P.cvholdout;
      end
      if Overwrite.lambda
        Queries.lambda = P.lambda;
      end
      if Overwrite.alpha
        Queries.alpha = P.alpha;
      end
      if ~isempty(Queries)
        qKeys = fieldnames(Queries);
        qVals = struct2cell(Queries);
        pKeys = fieldnames(P);
        pVals = struct2cell(P);
        z = ismember(pKeys, qKeys);
        pKeys = pKeys(z);
        pVals = pVals(z);
        pNum  = cellfun(@numel, pVals, 'Unif', 0);
        pNum  = cellfun(@(x) false(1,x), pNum, 'Unif', 0);
        Select = cell2struct(pNum, pKeys, 1);

        for q = 1:numel(qKeys)
          qk = qKeys{q};
          qv = qVals{q};
          if ~iscell(qv)
            qv = {qv};
          end
          try
            pv = P.(qk);
          catch ME
            fprintf('Error: Query key %s is not a field in Params.\n\n', qk)
            rethrow(ME);
          end
          if ~iscell(pv);
            pv = {pv};
          end
          if ischar(pv{1})
            if ~ischar(qv{1})
              error('Invalid type! The parameter field %s contains a string or cellstr.', qk);
            end
          elseif isnumeric(pv{1})
            if ~isnumeric(qv{1})
              error('Invalid type! The parameter field %s contains a number or numeric vector.', qk);
            end
          end
          % HACK
          if iscell(qv) && isnumeric(qv{1})
            qv = qv{1};
            pv = pv{1};
          end
          if any(strcmp(qk, {'data','cvholdout','lambda','alpha'}))
            Select.(qk) = ismember(pv, qv);
          else
            Select.(qk) = any(ismember(pv, qv));
          end
        end
        sKeys = fieldnames(Select);
        sVals = struct2cell(Select);
        DCLA = cell(1,4);
        dcla = {'data','cvholdout','lambda','alpha'};
        for k = 1:4
          DCLA{k} = Select.(dcla{k});
        end
        z = ~ismember({'data','cvholdout','lambda','alpha'}, sKeys);
        sVals = [DCLA,sVals(z)];
        sValExp = cell(1,numel(sVals));
        [sValExp{:}] = deal(ndgrid(sVals{:}));
        sValExp = cellfun(@(x) x(:), sValExp, 'Unif', 0);
        SEL = all(cell2mat(sValExp), 2);
      else
        SEL = true(N,1);
      end
    else
      SEL = true(N,1);
    end

    %% load results file
    resultpath = fullfile(jobdir, resultfile);
    if ~exist(resultpath, 'file') || ~any(SEL)
      continue;
    end
    tmp = load(resultpath);
    R = tmp.results;
    [R.job] = deal(i);

    for ii = 1:numel(R)
      iii = iii + 1;
      if SEL(ii)
        iiir = iiir + 1;
        tmp = R(ii);
        results(iiir) = tmp;
      end
    end
  end
  % clear out any empty entries.
  if iii < N;
    results((iii+1):N) = [];
  end
  fprintf('\n')
end

%% Private Functions
function y = selectcv(x,cv)
  if numel(x)==1
    if iscell(x);
      y = x{1};
    else
      y = x;
    end
    return
  end

  dim = size(x);
  if iscell(x)
    if length(dim)>1 && all(dim>1)
      y = x(cv,:);
    else
      y = x{cv};
    end
  else
    if length(dim)>1 && all(dim>1)
      y = x(cv,:);
    else
      y = x(cv);
    end
  end
end

function jobdirs = SelectJobDirs(dirs,varargin)
  p = inputParser;
  addRequired(p, 'dirs');
  addParameter(p, 'root','.', @ischar);
  addParameter(p, 'sort',false, @islogical);
  parse(p, dirs, varargin{:});

  dirs = p.Results.dirs;
  SORT = p.Results.sort;
  root = p.Results.root;

  N = length(dirs);
  isJobDir = false(N,1);
  for ii = 1:N
    jobdir = fullfile(root,dirs(ii).name);
    % Check if a special dir (current or parent)
    if any(strcmp(jobdir,{'.','..'}))
      continue
    end
    % Check if contains parameter file.
    paramfile = fullfile(jobdir, 'params.json');
    if exist(paramfile, 'file')
      isJobDir(ii) = true;
    end
  end
  jobdirs = dirs(isJobDir);

  if SORT
    try
      jobs = cellfun(@(x) sscanf(x,'%d'), {jobdirs.name});
      disp('Sorting jobs numerically.')
    catch
      jobs = {jobdirs.name};
      disp('Sorting jobs alphabetically.')
    end
    [~,ix] = sort(jobs);
    jobdirs = jobdirs(ix);
  end
  fprintf('Found %d job directories.\n', length(jobdirs))
end

function r = forceRowVec(v)
  r = v(1:end);
end

function c = forceColVec(v)
  c = v(:);
end
