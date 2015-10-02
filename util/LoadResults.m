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

  %% Loop over job dirs
  n = length(jobdirs);
  nchar = 0;
  fprintf('Loading job ');
  iii = 0;
  for i = 1:n;
    fprintf(repmat('\b', 1, nchar));
    nchar = fprintf('%d of %d', i, n);

    %% load parameter file
    jobdir      = fullfile(resultdir, jobdirs(i).name);
    parampath   = fullfile(jobdir,paramfile);
    tmp         = loadjson(parampath);
    tmp.jobdir  = jobdir;
    params(i)   = tmp;
    clear tmp;

    %% load results file
    resultpath = fullfile(jobdir, resultfile);
    if ~exist(resultpath, 'file')
      continue;
    end
    tmp = load(resultpath);
    R = tmp.results;
    [R.job] = deal(i);

    for ii = 1:numel(R)
      iii = iii + 1;
      results(iii) = R(ii);
    end
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
