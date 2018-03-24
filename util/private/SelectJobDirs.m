function jobDirs = SelectJobDirs(dirs, paramsFile, resultDir, sort)
  p = inputParser;
  addRequired(p, 'dirs');
  addRequired(p, 'ParamsFile', @ischar);
  addRequired(p, 'ResultDir', @ischar);
  addRequired(p, 'sort', @islogical);
  parse(p, dirs, paramsFile, resultDir, sort);

  JOB_DIRS = p.Results.dirs;
  PARAMS_FILE = p.Results.ParamsFile;
  RESULT_DIR = p.Results.ResultDir;
  SORT = p.Results.sort;

  N = length(JOB_DIRS);
  isJobDir = false(N,1);
  for ii = 1:N
    jobDir = fullfile(RESULT_DIR,JOB_DIRS(ii).name);
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
  jobDirs = {jobDirs.name};
  fprintf('Found %d job directories.\n', length(jobDirs))
end
