function [Params, Results, n] = HTCondorLoad(ResultDir, varargin)
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
    addParameter(p,'ResultFile','results.mat',@ischar);
    addParameter(p,'ParamFile','params.json',@ischar);
    addParameter(p,'SortJobs',false,@islogical)
    addParameter(p,'SkipFields',{})
    addParameter(p,'IncludeFields',{})
    addParameter(p,'JobList',{},@iscellstr)
    addParameter(p,'legacy',false,@islogical)
    addParameter(p,'quiet',false,@islogical)
    parse(p,ResultDir, varargin{:});

    RESULT_DIR  = p.Results.ResultDir;
    RESULT_FILE = p.Results.ResultFile;
    PARAMS_FILE = p.Results.ParamFile;
    SORT_JOBS   = p.Results.SortJobs;
    SKIP        = p.Results.SkipFields;
    INCLUDE     = p.Results.IncludeFields;
    QUIET       = p.Results.quiet;
    LEGACY      = p.Results.legacy;
    jobDirs     = p.Results.JobList;

    if isempty(jobDirs)
        jobDirs  = HTCondorListJobDirs(RESULT_DIR, PARAMS_FILE, 'SortJobs', SORT_JOBS);
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

%% Preallocate result structure
% Assumption is that all result files are identical size with same fields
    [Results,Params] = PreallocateResultsStruct(jobDirs,RESULT_DIR,RESULT_FILE,PARAMS_FILE);
    if isfield(Results,'nz_rows') && ~isfield(Results,'Uix')
        Results(1).Uix = [];
    end
    Results(1).jobdir = [];
    Params(1).subject = [];
    Params(1).jobdir = [];
    nchar = 0;
    if ~QUIET
        fprintf('Loading job ');
    end
    cursor = 0;
    n = zeros(nJobDirs,1);
    for i = 1:nJobDirs;
        if ~QUIET
            fprintf(repmat('\b', 1, nchar));
            nchar = fprintf('%d of %d', i, nJobDirs);
        end

        % load parameter file
        jobDir      = fullfile(RESULT_DIR, jobDirs{i});
        paramsPath  = fullfile(jobDir,PARAMS_FILE);
        tmp         = loadjson(paramsPath);
        if iscell(tmp.data)
            tmp.subject = cellfun(@(x) sscanf(x,'s%02d'), tmp.data, 'Unif', 0);
        else
            tmp.subject = sscanf(tmp.data,'s%02d');
        end
        tmp.jobdir  = jobDir;
        tmp = orderfields(tmp, Params(1));
        Params(i) = tmp;
        clear tmp;

        if nargout > 1
            % load results file
            resultPath = fullfile(jobDir, RESULT_FILE);
            if ~exist(resultPath, 'file')
                continue;
            end
            tmp = load(resultPath);
            R = tmp.results;
            [R.jobdir] = deal(i);
            R = rmfield(R, SKIP);
            a = cursor + 1;
            b = cursor + numel(R);

            if isfield(R,'nz_rows') && ~isfield(R,'Uix')
                for iii = 1:numel(R)
                    R(iii).Uix = find(R(iii).nz_rows);
                end
%                 z = ismember(fieldnames(R),fieldnames(Results));
%                 fnp = fieldnames(R);
%                 fnp(z) = [];
                R = orderfields(R, Results(1));
                Results(a:b) = R;

                
                
            elseif isfield(R,'Uz')
                for iii = 1:numel(R)
                    R(iii).Sz = [];
                    R(iii).nz_rows = full(any(R(iii).Uz, 2));
                    R(iii).Uix = find(R(iii).nz_rows);
                    R(iii).nvox = size(R(iii).Uz, 1);
                    R(iii).structureScoreMap = [];
                end
                if a == 1
                    [Results.Sz] = deal([]);
                    [Results.nz_rows] = deal([]);
                    [Results.Uix] = deal([]);
                    [Results.nvox] = deal([]);
                end
                z = ismember(fieldnames(R),fieldnames(Results));
                fnp = fieldnames(R);
                fnp(z) = [];
                R = orderfields(R, rmfield(R,fnp{:}));
                Results(a:b) = R;
                %               Results(a:b) = R;
            else
                z = ismember(fieldnames(R),fieldnames(Results));
                fnp = fieldnames(R);
                fnp(z) = [];
                if all(z)
                    Results(a:b) = orderfields(R,Results(1));
                else
                    Results(a:b) = orderfields(rmfield(R,fnp{:}),Results(1));
                end
            end
            n(i) = numel(R);
            cursor = b;
        end
    end
    
    if ~QUIET
        fprintf('\n')
    end
end

function [results,params] = PreallocateResultsStruct(jobDirs,RESULT_DIR,RESULT_FILE,PARAMS_FILE)
    [r,p] = LoadFirstNotEmpty(jobDirs,RESULT_DIR,RESULT_FILE,PARAMS_FILE);
    results = repmat(EmptyFields(r(:)),numel(jobDirs),1);
    params = repmat(EmptyFields(p(:)),numel(jobDirs),1);
end

function [r,p] = LoadFirstNotEmpty(jobDirs,RESULT_DIR,RESULT_FILE,PARAMS_FILE)
    i = 0;
    fileExists = false;
    while ~fileExists
        i = i + 1;
        jobDir     = fullfile(RESULT_DIR, jobDirs{i});
        resultPath = fullfile(jobDir, RESULT_FILE);
        if exist(resultPath, 'file') == 2
            tmp = load(resultPath);
            fileExists = isfield(tmp,'results');
        end
    end
    r = tmp.results;
    if nargout > 1
        paramsPath = fullfile(jobDir, PARAMS_FILE);
        p = loadjson(paramsPath);
    end
end

function [se] = EmptyFields(s)
    f = fieldnames(s);
    x = [f(:)';repmat({[]},1,numel(f))];
    se = struct(x{:});
end
    

% if ~isfield(tmp,'results')
%     tmp2 = tmp; clear tmp;
%     m = numel(tmp2.err1);
%     tmp.results = init_results();
%     tmp.results(m).job = [];
%     clear tmp2;
% end