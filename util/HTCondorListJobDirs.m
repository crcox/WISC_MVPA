function [ jobDirs ] = HTCondorListJobDirs( varargin )
    p = inputParser();
    addOptional(p,'ResultDir', pwd(), @ischar);
    addOptional(p,'ParamFile','params.json',@ischar);
    addParameter(p,'SortJobs',false,@islogical)
    parse(p,varargin{:});

    basedir = GetFullPath(p.Results.ResultDir);
    if ~exist(basedir,'dir')
        error('There does not seem to be a directory at %s ...', basedir);
    end
    PARAMS_FILE = p.Results.ParamFile;
    SORT_JOBS   = p.Results.SortJobs;
    
    allFiles = dir(basedir);
    allDirs  = allFiles([allFiles.isdir]);
    jobDirs  = SelectJobDirs(allDirs, PARAMS_FILE, basedir, SORT_JOBS);
end

