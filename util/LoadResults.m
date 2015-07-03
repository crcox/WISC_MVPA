function [results, params] = LoadResults(varargin)
  p = inputParser();
  addParameter(p,'ResultDir','.',@ischar);
  addParameter(p,'DataDir','/media/chris/FacePlaceObject/SOSLASSO/SOSLasso/batch-002/tuning', @ischar)
  addParameter(p,'MetadataFile','metadata.mat',@ischar)
  addParameter(p,'Filters',{'rowfilter'},@iscell)
  addParameter(p,'SortJobs',false,@islogical)
  addParameter(p,'SkipFields',[])
  addParameter(p,'ResultFile','fitObj.mat',@ischar);
  addParameter(p,'ParamFile','params.json',@ischar);
  parse(p,varargin{:});

  resultdir  = p.Results.ResultDir;
  datadir    = p.Results.DataDir;
  metafile   = p.Results.MetadataFile;
  filters    = p.Results.Filters;
  sortjobs   = p.Results.SortJobs;
  SKIP       = p.Results.SkipFields;
  resultfile = p.Results.ResultFile;
  paramfile  = p.Results.ParamFile;
  allfiles   = dir(resultdir);
  alldirs    = allfiles([allfiles.isdir]);
  jobdirs    = SelectJobDirs(alldirs, paramfile, sortjobs);

  if ~isempty(SKIP)
    SkipStr = SKIP;
    SkipStr{end} = ['and ', SKIP{end}];
    SkipStr = strjoin(SkipStr, ', ');
    fprintf('Fields %s will be skipped.\n', SkipStr);
  end

  metapath = fullfile(datadir,metafile);
  load(metapath, 'metadata');
  cvscheme = 1;

  n = length(jobdirs);
  nchar = 0;
  fprintf('Loading job ');
  ind = 0;
  for i = 1:n;
    fprintf(repmat('\b', 1, nchar));
    nchar = fprintf('%d of %d', i, n);

    jobdir      = jobdirs(i).name;
    parampath   = fullfile(jobdir,paramfile);
    tmp         = loadjson(parampath);
    tmp.jobdir  = jobdir;
    params(i)   = tmp;

    %subjectNumbers = cellfun(@(x) sscanf(x, 'jlp%02d.mat'), params.data);
    subjectNumbers = 1:10;
    includedSubjects = ismember([metadata.subject], subjectNumbers);
    metadata = metadata(includedSubjects);
    [~,sidx] = ismember(subjectNumbers, [metadata.subject]);

    cvIndexes = 1:size(metadata(1).CVBLOCKS,2);
    finalind = params.FinalHoldoutSet;
    cvIndexes(finalind) = [];

    % These will be used when looping over the cv and subject elements of the
    % fitObj.
    [snum, cvind] = ndgrid(subjectNumbers,cvIndexes);

    filter = cell(length(subjectNumbers),1);
    cvfilter = cell(length(subjectNumbers),length(cvIndexes));
    finalfilter = cell(length(subjectNumbers),length(cvIndexes));
    for ss = 1:length(subjectNumbers);
      for ii = 1:length(filters);
        fname = filters{ii};
        z = metadata(sidx(ss)).(fname);
        if ii == 1;
          filter{ss} = z;
        else
          filter{ss} = filter{ss} & z;
        end
      end
      for cc = 1:length(cvIndexes)
        cvfilter{ss,cc} = metadata(ss).CVBLOCKS(filter{ss},cvIndexes(cc));
        finalfilter{ss,cc} = metadata(ss).CVBLOCKS(filter{ss},finalind);
      end
    end
%    cvind       = params(i).cvholdout;
%    finalind    = params(i).finalholdout;
%    cvscheme    = params(i).cvscheme;

    resultpath = fullfile(jobdir, resultfile);
    if exist(resultpath, 'file')
      load(resultpath, 'fitObj');
      N = numel(fitObj);
      for ii = 1:N
        train = ~cvfilter{ii}(~finalfilter{ii});
        test = cvfilter{ii}(~finalfilter{ii});
        ind = ind + 1;
        tmp             = fitObj(ii);
        tmp.dp1         = computeStat(fitObj(ii), test, 'dprime');
        tmp.dp2         = computeStat(fitObj(ii), train, 'dprime');
        tmp.err1        = computeStat(fitObj(ii), test, 'error');
        tmp.err2        = computeStat(fitObj(ii), train, 'error');
        tmp.nz_rows     = computeStat(fitObj(ii), [], 'nzvox');
        tmp.cvfilter    = cvfilter;
        tmp.finalfilter = finalfilter;
        tmp.pind        = i;
        tmp.cvholdout   = cvind(ii);
        tmp.subject     = snum(ii);
        if ~isempty(SKIP)
          nskip = length(SKIP);
          for iii = 1:nskip
            tmp = rmfield(tmp, SKIP{iii});
          end
        end
        results(ind)  = tmp;
      end
    end
  end
  fprintf(' %d\n', ind)
end

function x = computeStat(fitObj, filter, stat)
  y  = fitObj.Y(filter) > 0;
  yz = fitObj.Yh(filter) > 0;
  switch stat
  case 'error'
    x = 1-(sum(y==yz)/length(y));
  case 'dprime'
    ntp = sum(y);
    nfp = sum(~y);
    x = (sum(y&yz)/ntp) - (sum(~y&yz)/nfp);
  case 'nzvox'
    x = sum(fitObj.betas~=0);
  end
end

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
function jobdirs = SelectJobDirs(dirs,paramfile, SORT)
  N = length(dirs);
  isJobDir = false(N,1);
  for ii = 1:N
    jobdir = dirs(ii).name;
    % Check if a special dir (current or parent)
    if any(strcmp(jobdir,{'.','..'}))
      continue
    end
    % Check if contains parameter file.
    parampath = fullfile(jobdir, paramfile);
    if exist(parampath, 'file')
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
