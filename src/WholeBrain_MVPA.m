function WholeBrain_MVPA(varargin)
  p = inputParser;
  p.KeepUnmatched = false;
  %% Parse and set parameters
  addParameter(p , 'debug'            , false     , @islogicallike );
  addParameter(p , 'SmallFootprint'   , false     , @islogicallike );
  addParameter(p , 'regularization'        , []        , @ischar        );
  addParameter(p , 'debias'           , false      , @islogicallike );
  addParameter(p , 'normalize'        , false                      );
  addParameter(p , 'bias'             , false     , @islogicallike );
  addParameter(p , 'filters'          , []        , @ischarlike    );
  addParameter(p , 'target'           , []        , @ischar        );
  addParameter(p , 'data'             , []        , @ischarlike    );
  addParameter(p , 'data_var'         , 'X'       , @ischar        );
  addParameter(p , 'metadata'         , []        , @ischar        );
  addParameter(p , 'metadata_var'     , 'metadata', @ischar        );
  addParameter(p , 'cvscheme'         , []        , @isintegerlike );
  addParameter(p , 'cvholdout'        , []        , @isnumeric     );
  addParameter(p , 'orientation'      , []        , @ischar        );
  addParameter(p , 'diameter'         , []        , @isnumeric     );
  addParameter(p , 'overlap'          , []        , @isnumeric     );
  addParameter(p , 'shape'            , []        , @ischar        );
  addParameter(p , 'slradius'         , []        , @isnumeric     );
  addParameter(p , 'slTestToUse'      , 'accuracyOneSided_analytical', @ischar);
  addParameter(p , 'slpermutations'   , 0         , @isscalar      );
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
  addParameter(p , 'RandomSeed'       , 0                          );
  addParameter(p , 'PermutationTest'  , false     , @islogicallike );
  addParameter(p , 'SaveResultsAs'  , 'mat'     , @isMatOrJSON);

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
  required = {'regularization','data','metadata','cvscheme','cvholdout','finalholdout'};
  assertRequiredParameters(p.Results,required);

  DEBUG            = p.Results.debug;
  SmallFootprint   = p.Results.SmallFootprint;
  regularization        = p.Results.regularization;
  debias           = p.Results.debias;
  normalize        = p.Results.normalize;
  BIAS             = p.Results.bias;
  filter_labels    = p.Results.filters;
  target_label     = p.Results.target;
  datafile         = p.Results.data;
  data_var         = p.Results.data_var;
  cvscheme         = p.Results.cvscheme;
  cvholdout        = p.Results.cvholdout;
  orientation      = p.Results.orientation;
  diameter         = p.Results.diameter;
  overlap          = p.Results.overlap;
  shape            = p.Results.shape;
  slradius         = p.Results.slradius;
  slTestToUse      = p.Results.slTestToUse;
  slpermutations   = p.Results.slpermutations;
  finalholdoutInd  = p.Results.finalholdout;
  metafile         = p.Results.metadata;
  metadata_var     = p.Results.metadata_var;
  lambda           = p.Results.lambda;
  alpha            = p.Results.alpha;
  opts             = p.Results.AdlasOpts;
  SanityCheckData  = p.Results.SanityCheckData; %#ok<NASGU>
  RandomSeed       = p.Results.RandomSeed;
  PermutationTest  = p.Results.PermutationTest;
  SaveResultsAs    = p.Results.SaveResultsAs;

  rng(RandomSeed);

  p.Results

  % Check that the correct parameters are passed, given the desired regularization
  [lambda, alpha] = verifyLambdaSetup(regularization, lambda, alpha);

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
      [rowfilter{i},colfilter{i}] = composeFilters(metadata(i).filters, filter_labels);
      if isempty(rowfilter{i})
        rowfilter{i} = true(1,metadata(i).nrow);
      end
      if isempty(colfilter{i})
        colfilter{i} = true(1,metadata(i).ncol);
      end
    end
  end

  %% Load CV indexes, and identify the final holdout set.
  % N.B. the final holdout set is excluded from the rowfilter.
  cvind = cell(1,N);
  cvindAll = cell(1,N);
  for i = 1:N
    % Add the final holdout set to the rowfilter, so we don't even load
    % those data.
    cvindAll{i} = metadata(i).cvind(:,cvscheme);
    finalholdout = cvindAll{i} == finalholdoutInd;
    % Remove the final holdout set from the cvind, to match.
    rowfilter{i} = forceRowVec(rowfilter{i}) & forceRowVec(~finalholdout);
    cvind{i} = cvindAll{i}(rowfilter{i});
  end

  %% Load data and select targets
  [X,subjix] = loadData(datafile, data_var, rowfilter, colfilter, metadata);
  metadata   = metadata(subjix);
  rowfilter  = rowfilter(subjix);
  colfilter  = colfilter(subjix);
  cvind      = cvind(subjix);

  %% Select targets
  fprintf('\n');
  fprintf('Loading similarity structure\n');
  fprintf('----------------------------\n');
  fprintf('%12s: %s\n', 'target_label', target_label);
  fprintf('%12s: %s\n', 'type', 'category');
  fprintf('\n');
  Y = selectTargets(metadata, 'category', target_label, [], [], rowfilter);

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
    msg = 'YES';
    X = ascell(X);
    for ii = 1:numel(X)
      switch normalize
        case 'zscore'
          mm = mean(X{ii},1);
          ss = std(X{ii},0,1);
        case 'stdev'
          mm = 0;
          ss = std(X{ii},0,1);
        case '2norm'
          mm = mean(X{ii},1);
          ss = norm(X{ii});
        otherwise
          % N.B. If zscore_train, that is handled later.
          mm = 0;
          ss = 1;
      end
      X{ii} = bsxfun(@minus,X{ii}, mm);
      X{ii} = bsxfun(@rdivide,X{ii}, ss);
    end
  end
  fprintf(': [%3s]\n', msg);

  % Final holdout index
  fprintf('%-26s', 'Final holdout index');
  fprintf(': [%3d]\n', finalholdoutInd);

  fprintf('Data loaded and processed.\n');

  %% Plug in the parameters and run
  switch regularization
  case 'lasso'
    [results,info] = learn_category_encoding(Y, X, regularization, ...
                      'lambda'         , lambda         , ...
                      'alpha'          , alpha          , ...
                      'cvind'          , cvind          , ...
                      'cvholdout'      , cvholdout      , ...
                      'normalize'      , normalize      , ...
                      'DEBUG'          , DEBUG          , ...
                      'SmallFootprint' , SmallFootprint , ...
                      'PermutationTest', PermutationTest, ...
                      'AdlasOpts'      , opts); %#ok<ASGLU>

  case 'searchlight'
    X = uncell(X);
    Y = uncell(Y)+1;
    cvind = uncell(cvind);
    colfilter = uncell(colfilter);

    % create a 3D binary mask
    z = strcmp({metadata.coords.orientation}, orientation);
    xyz = metadata.coords(z).xyz(colfilter,:);
    [mask,dxyz] = coordsTo3dMask(xyz);

    % Translate slradius (in mm) to sl voxels
    % N.B. Because voxels need not be symmetric cubes, but Seachmight will
    % generate symmetric spheres from a single radius parameter, we need to
    % select one value of the three that will be produced in this step. I am
    % arbitrarily choosing the max, to err on the side of being inclusive.
    slradius_ijk = max(round(slradius ./ dxyz));

    % create the "meta" neighbourhood structure
    meta = createMetaFromMask(mask, slradius_ijk);

    % Prepare parameters
    classifier = 'gnb_searchmight';
    if strcmp(slTestToUse,'accuracyOneSided_permutation')
      TestToUseCfg = {'testToUse',slTestToUse,slpermutations};
    else
      TestToUseCfg = {'testToUse',slTestToUse};
    end
    [am,pm,hm,fm] = computeInformationMap(X,Y,cvind,classifier,'searchlight', ...
                                meta.voxelsToNeighbours,meta.numberOfNeighbours,TestToUseCfg{:});

    results.accuracy_map = am;
    results.hitrate_map = hm;
    results.falsealarm_map = fm;
    results.pvalue_map = pm;

  case 'soslasso'
    xyz = cell(numel(metadata),1);
    for ii = 1:numel(xyz)
      z = strcmp({metadata(ii).coords.orientation}, orientation);
      xyz{ii} = metadata(ii).coords(z).xyz(colfilter{ii},:);
    end
    G = coordGrouping(xyz, diameter, overlap, shape);
    [results,info] = learn_category_encoding(Y, X, regularization, ...
                      'groups'         , G              , ...
                      'lambda'         , lambda         , ...
                      'alpha'          , alpha          , ...
                      'cvind'          , cvind          , ...
                      'cvholdout'      , cvholdout      , ...
                      'normalize'      , normalize      , ...
                      'DEBUG'          , DEBUG          , ...
                      'debias'         , debias         , ...
                      'SmallFootprint' , SmallFootprint , ...
                      'PermutationTest', PermutationTest, ...
                      'AdlasOpts'      , opts); %#ok<ASGLU>
    %% Revise cv indexes
    % Add the final holdout index to all results.
    [results.finalholdout] = deal(finalholdoutInd);
    % Adjust the cvholdout indexes to accomodate the final holdout index.
    if isfield(results,'cvholdout') && finalholdoutInd > 0
      cvholdout = [results.cvholdout];
      z = cvholdout >= finalholdoutInd;
      cvholdout(z) = cvholdout(z) + 1;
      cvholdout = mat2cell(cvholdout(:),ones(numel(cvholdout),1));
      [results.cvholdout] = deal(cvholdout{:});
    end
    %% Add extra parameter info
    [results.diameter] = deal(diameter);
    [results.overlap] = deal(overlap);
    [results.shape] = deal(shape);

  end
  whos results
  fprintf('Saving %d results\n', numel(results));
  fprintf('\t%s\n', 'results.mat');

  %% Save results
  rinfo = whos('results');
  switch SaveResultsAs
      case 'mat'
          if rinfo.bytes > 2e+9
            save('results.mat','results','-v7.3');
          else
            save('results.mat','results');
          end
      case 'json'
          savejson('',results,'FileName','results.json','ForceRootName',false);
  end
  fprintf('Done!\n');
end

%% Local functions
function [lambda, alpha] = verifyLambdaSetup(regularization, lambda, alpha)
% Each regularization requires different lambda configurations. This private
% function ensures that everything has been properly specified.
  switch regularization
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

function b = isMatOrJSON(x)
    b = any(strcmpi(x, {'mat','json'}));
end
