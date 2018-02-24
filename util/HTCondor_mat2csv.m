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
    addParameter(p,'correct_nzv',0,@isnumeric)
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
    correct_nzv = p.Results.correct_nzv;
    QUIET       = p.Results.quiet;
    jobDirs     = p.Results.JobList;

    if isempty(jobDirs)
        allFiles = dir(RESULT_DIR);
        allDirs  = allFiles([allFiles.isdir]);
        jobDirs  = SelectJobDirs(fullfile(RESULT_DIR,{allDirs.name}), PARAMS_FILE, SORT_JOBS);
    else
        jobDirs = fullfile(RESULT_DIR, jobDirs);
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
    z = ismember(SKIP,fieldnames(R));
    R = rmfield(R, SKIP(z));
    [fmt,fmt_h] = fieldfmt(fieldnames(R));

    fid = fopen(OUTPUT_FILE, 'w');
    header = fieldnames(R);
    fprintf(fid, fmt_h, header{:});
    nchar = 0;
    SHOW_WARNING = 1;
    for i = 1:nJobDirs
        if ~QUIET
            fprintf(repmat('\b', 1, nchar));
            nchar = fprintf('%d of %d', i, nJobDirs);
        end

        % load parameter file
        jobDir = jobDirs{i}; %fullfile(RESULT_DIR, jobDirs{i});
        pfile  = fullfile(jobDir,PARAMS_FILE);
        rfile  = fullfile(jobDir,RESULT_FILE);
        PARAMS_LOADED = 0;
        if exist(rfile, 'file') && exist(pfile, 'file')
            r = load(rfile, 'results');
            if isempty(r.results(1).subject)
                p = loadjson(pfile);
                PARAMS_LOADED = 1;
                if iscell(p.data) && numel(p.data) > 1 && SHOW_WARNING
                    warning('Subject ID is missing from results structure and params specifies multiple data files. Cannot auto-populate subject IDs');
                    SHOW_WARNING = 0;
                else
                    subjectID = extractSubjectID(p.data, p.subject_id_fmt);
                end
            end
            z = ismember(SKIP,fieldnames(r.results));
            R = rmfield(r.results, SKIP(z));
            if isempty(r.results(1).subject)
                [R.subject] = deal(subjectID);
            end
            if correct_nzv
                for ii = 1:numel(R)
                    l2norm = sum(r.results(ii).Uz .^ 2, 2);
                    R(ii).nzv = nnz(l2norm > correct_nzv);
                end
            end
            if isfield(R,'iterations')
                for ii = 1:numel(R)
                    R(ii).iterations = numel(R(ii).iterations);
                end
            end
            if isfield(R,'RandomSeed')
                for ii = 1:numel(R)
                    if (numel(R(ii).RandomSeed) > numel(R(ii).err1)) && (numel(R(ii).RandomSeed) == numel(R))
                        R(ii).RandomSeed = R(ii).RandomSeed(ii);
                    end
                end
            end
            % BeyondMagnitude Hack
            if isfield(R,'subject')
                for ii = 1:numel(R)
                    if ischar(R(ii).subject);
                        R(ii).subject = sscanf(R(ii).subject, 'BM%d.mat');
                    end
                end
            end
            if isfield(R,'nzv')&& isempty(R(1).nzv) && isfield(R,'nz_rows') && ~isempty(R(1).nz_rows)
                for ii = 1:numel(R)
                    R(ii).nzv = nnz(R(ii).nz_rows);
                    R(ii).nvox = numel(R(ii).nz_rows);
                end
            end
            if isfield(R,'finalholdout')&& islogical(R(1).finalholdout)
                if ~PARAMS_LOADED
                    p = loadjson(pfile);
%                     PARAMS_LOADED = 1;
                end
                for ii = 1:numel(R)
                    R(ii).finalholdout = p.finalholdout;
                end  
            end
            rc = squeeze(struct2cell(R));
            
            for j = 1:size(rc,2)
                xc = rc(:,j);
                for k = 1:size(rc,1);
                    if isinteger(rc(k,j))
                        xc{k} = double(rc(k,j));
                    end
                end
%                 fprintf(fmt, xc{:})
                fprintf(fid, fmt, xc{:});
            end
        end
    end
%     fprintf('\n');
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
        'cvscheme'     , '%d'   , ...
        'data'         , '%s'   , ...
        'data_var'     , '%s'   , ...
        'data_varname' , '%s'   , ...
        'metadata'     , '%s'   , ...
        'metadata_var' , '%s'   , ...
        'metadata_varname' , '%s'   , ...
        'target'       , '%s'   , ...
        'target_type'  , '%s'   , ...
        'target_label' , '%s'   , ...
        'sim_source'   , '%s'   , ...
        'sim_metric'   , '%s'   , ...
        'Gtype'        , '%s'   , ...
        'regularization', '%s'   , ...
        'shape'        , '%s'   , ...
        'alpha'        , '%.4f' , ...
        'lambda'       , '%.4f' , ...
        'lambda1'      , '%.4f' , ...
        'LambdaSeq'    , '%s'   , ...
        'diameter'     , '%d'   , ...
        'overlap'     , '%d'   , ...
        'tau'          , '%.4f' , ...
        'normalize'    , '%s'   , ...
        'normalize_data', '%s'   , ...
        'normalize_target', '%s'   , ...
        'normalizewrt', '%s'   , ...
        'target_normalization'    , '%s'   , ...
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
        'iterations'   , '%d' , ...
        'p1'           , '%.4f' , ...
        'p2'           , '%.4f' , ...
        'cor1'         , '%.4f' , ...
        'cor2'         , '%.4f' , ...
        'err1'         , '%.4f' , ...
        'err2'         , '%.4f' , ...
        'FroErr1'      , '%.4f' , ...
        'FroErr2'      , '%.4f' , ...
        'nvox'         , '%d'   , ...
        'Wnz'          , '%d'   , ...
        'nzv'          , '%d');

    fmt_c = cell(1,numel(fieldnames));
    for i = 1:numel(fieldnames)
        fmt_c{i} = fieldFMT.(fieldnames{i});
    end
    fmt = [strjoin(fmt_c,',') '\n'];
    fmt_h = [strjoin(repmat({'%s'},1,numel(fieldnames)),',') '\n'];
end
