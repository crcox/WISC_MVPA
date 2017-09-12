function HTCondor_mat2csv(ResultDir, varargin)
% HTCondorLoad Load results from a collection of HTCondor jobs.
%   [params] = HTCondorLoad(ResultDir); avoids loading any data.
%   [params, results] = HTCondorLoad(ResultDir);
%   [params, results] = HTCondorLoad(ResultDir, varargin);
%
%   OPTIONS:
%   ResultFile : The name of the file written by the HTCondor job. Will
%                contain a structure.
%   ParamFile  : The name of the json file that defines the job.
%   SortJobs   : Will ensure that the internal list of job directories is
%                sorted by directory name before loading any data.
%   SkipFields : A cell array of field names in the results structure. Any
%                field name listed here will not be loaded. Useful if there
%                are large data elements that you do not need. Controls
%                memory footprint and speeds access.
%   JobList    : A list of job directory names to load.
%
%  DEPENDENCIES:
%  - JSONlab http://www.mathworks.com/matlabcentral/fileexchange/33381
    p = inputParser();
    addRequired(p,'ResultDir', @ischar);
    addRequired(p,'OutputFile', @ischar);
    addParameter(p,'ResultFile','results.mat',@ischar);
    addParameter(p,'ParamFile','params.json',@ischar);
    addParameter(p,'SortJobs',false,@islogical)
    addParameter(p,'SkipFields',{})
    addParameter(p,'IncludeFields',{})
    addParameter(p,'JobList',{},@iscellstr)
    addParameter(p,'quiet',false,@islogical)
    parse(p,ResultDir, varargin{:});

    RESULT_DIR  = p.Results.ResultDir;
    RESULT_FILE = p.Results.ResultFile;
    PARAMS_FILE = p.Results.ParamFile;
    OUTPUT_FILE = p.Results.OutputFile;
    SORT_JOBS   = p.Results.SortJobs;
    SKIP        = p.Results.SkipFields;
    INCLUDE     = p.Results.IncludeFields;
    QUIET       = p.Results.quiet;
    jobDirs     = p.Results.JobList;

    if isempty(jobDirs)
        allFiles = dir(RESULT_DIR);
        allDirs  = allFiles([allFiles.isdir]);
        jobDirs  = SelectJobDirs(fullfile(RESULT_DIR,{allDirs.name}), PARAMS_FILE, SORT_JOBS);
    end
    nJobDirs = numel(jobDirs);

    if ~isempty(SKIP) && ~isempty(INCLUDE)
        error('Both SKIP and INCLUDE cannot be set.')
    elseif ~isempty(SKIP) && isempty(INCLUDE)
        SkipStr = SKIP;
        SkipStr{end} = ['and ', SKIP{end}];
        SkipStr = strjoin(SkipStr, ', ');
        fprintf('Fields %s will be skipped.\n', SkipStr);
    elseif isempty(SKIP) && ~isempty(INCLUDE)
        InclStr = INCLUDE;
        InclStr{end} = ['and ', INCLUDE{end}];
        InclStr = strjoin(InclStr, ', ');
        fprintf('Fields %s will be Included.\n', InclStr);
    end

    i = 0;
    fileExists = false;
    while ~fileExists
        i = i + 1;
        jobDir     = jobDirs{i}; %fullfile(RESULT_DIR, jobDirs{i});
        rfile = fullfile(jobDir, RESULT_FILE);
        fileExists = exist(rfile, 'file') == 2;
    end
    tmp = load(rfile);
    R = tmp.results;
    if ~isempty(INCLUDE) && isempty(SKIP)
        SKIP = fieldnames(rmfield(R,INCLUDE));
    end
    R = rmfield(R, SKIP);
    [fmt,fmt_h] = fieldfmt(fieldnames(R));

    fid = fopen(OUTPUT_FILE, 'w');
    header = fieldnames(R);
    fprintf(fid, fmt_h, header{:});
    nchar = 0;
    for i = 1:nJobDirs
        if ~QUIET
            fprintf(repmat('\b', 1, nchar));
            nchar = fprintf('%d of %d', i, nJobDirs);
        end

        % load parameter file
        jobDir = jobDirs{i}; %fullfile(RESULT_DIR, jobDirs{i});
        pfile  = fullfile(jobDir,PARAMS_FILE);
        rfile  = fullfile(jobDir,RESULT_FILE);
        if exist(rfile, 'file') && exist(pfile, 'file')
            r = load(rfile, 'results');
            R = rmfield(r.results, SKIP);
            rc = squeeze(struct2cell(R));
            
            for j = 1:size(rc,2)
                xc = rc(:,j);
                for k = 1:size(rc,1);
                    if isinteger(rc(k,j))
                        xc{k} = double(rc(k,j));
                    end
                end
                fprintf(fid, fmt, xc{:});
            end
        end
    end
    fprintf('\n');
    fclose(fid);
end
function jobDirs = SelectJobDirs(dirs, paramsFile, sort)
  p = inputParser;
  addRequired(p, 'dirs');
  addRequired(p, 'ParamsFile', @ischar);
  addRequired(p, 'sort', @islogical);
  parse(p, dirs, paramsFile, sort);

  JOB_DIRS = p.Results.dirs;
  PARAMS_FILE = p.Results.ParamsFile;
  SORT = p.Results.sort;

  N = length(JOB_DIRS);
  isJobDir = false(N,1);
  for ii = 1:N
    jobDir = JOB_DIRS{ii}; %(ii).name;
    % Check if a special dir (current or parent)
    if any(strcmp(jobDir,{'.','..'}))
      continue
    end
    % Check if contains parameter file.
    paramsPath = fullfile(jobDir, PARAMS_FILE);
    if exist(paramsPath, 'file')
      isJobDir(ii) = true;
    end
  end
  jobDirs = JOB_DIRS(isJobDir);

  if SORT
    try
      jobs = cellfun(@(x) sscanf(x,'%d'), {jobDirs.name});
      disp('Sorting jobs numerically.')
    catch
      jobs = {jobDirs.name};
      disp('Sorting jobs alphabetically.')
    end
    [~,ix] = sort(jobs);
    jobDirs = jobDirs(ix);
  end
%   jobDirs = {jobDirs.name};
  fprintf('Found %d job directories.\n', length(jobDirs))
end

function [fmt,fmt_h] = fieldfmt(fieldnames)
    % List all possible fields and their desired format code
    fieldFMT = struct( ...
        'iter'         , '%d'   , ...
        'job'          , '%d'   , ...
        'subject'      , '%d'   , ...
        'finalholdout' , '%d'   , ...
        'cvholdout'    , '%d'   , ...
        'data'         , '%s'   , ...
        'target'       , '%s'   , ...
        'Gtype'        , '%s'   , ...
        'regularization', '%s'   , ...
        'lambda'       , '%.4f' , ...
        'lambda1'      , '%.4f' , ...
        'LambdaSeq'    , '%s'   , ...
        'tau'          , '%.4f' , ...
        'normalize'    , '%s'   , ...
        'bias'         , '%d'   , ...
        'RandomSeed'   , '%d'   , ...
        'nVoxel'       , '%d'   , ...
        'h1'           , '%d' , ...
        'h2'           , '%d' , ...
        'f1'           , '%d' , ...
        'f2'           , '%d' , ...
        'nt1'          , '%d' , ...
        'nt2'          , '%d' , ...
        'nd1'          , '%d' , ...
        'nd2'          , '%d' , ...
        'cor1'         , '%.4f' , ...
        'cor2'         , '%.4f' , ...
        'err1'         , '%.4f' , ...
        'err2'         , '%.4f' , ...
        'FroErr1'      , '%.4f' , ...
        'FroErr2'      , '%.4f' , ...
        'nvox'         , '%d'   , ...
        'nzv'          , '%d');

    fmt_c = cell(1,numel(fieldnames));
    for i = 1:numel(fieldnames)
        fmt_c{i} = fieldFMT.(fieldnames{i});
    end
    fmt = [strjoin(fmt_c,',') '\n'];
    fmt_h = [strjoin(repmat({'%s'},1,numel(fieldnames)),',') '\n'];
end
