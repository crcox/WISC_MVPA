function WholeBrain_RSA(varargin)
    p = inputParser;
    p.KeepUnmatched = false;
    % ----------------------Set parameters-----------------------------------------------
    % Model definition
	addParameter(p , 'regularization' , []        , @ischar        );
    addParameter(p , 'bias'             , false   , @islogicallike );
    addParameter(p , 'alpha'            , []      , @isnumeric     );
    addParameter(p , 'lambda'           , []      , @isnumeric     );
    addParameter(p , 'lambda1'          , []      , @isnumeric     );
    addParameter(p , 'lamSOS'           , []      , @isnumeric     );
    addParameter(p , 'lamL1'            , []      , @isnumeric     );
    addParameter(p , 'LambdaSeq'        , []      , @ischar        );
    addParameter(p , 'AdlasOpts'        , struct(), @isstruct      );
    addParameter(p , 'diameter'         , []      , @isnumeric     );
    addParameter(p , 'shape'            , []      , @ischar        );
    addParameter(p , 'overlap'          , []      , @isnumeric     );
    % Target definition
    addParameter(p , 'target_label'     , [], @ischar );
    addParameter(p , 'target_type'      , [], @ischar );
    addParameter(p , 'sim_source'       , [], @ischar );
    addParameter(p , 'sim_metric'       , [], @ischar );
    addParameter(p , 'tau'              , [], @isnumeric     );
    % Data definition
    addParameter(p , 'filters'          , []                  );
    addParameter(p , 'data'             , []                  );
    addParameter(p , 'data_varname'     , []                  );
    addParameter(p , 'metadata'         , [] , @ischar        );
    addParameter(p , 'metadata_varname' , [] , @ischar        );
    addParameter(p , 'finalholdout'     , 0  , @isintegerlike );
    addParameter(p , 'cvscheme'         , [] , @isnumeric     );
    addParameter(p , 'cvholdout'        , [] , @isnumeric     );
    addParameter(p , 'orientation'      , [] , @ischar        );
    % Normalization
    addParameter(p , 'normalize_data'  , 'none'        , @ischar ); % NEW NAME!
    addParameter(p , 'normalize_target', 'none'        , @ischar ); % NEW VARIABLE!
    addParameter(p , 'normalize_wrt'    , 'all_examples', @ischar ); % NEW VARIABLE!
    % Permutation
    addParameter(p , 'RandomSeed'             , []                    );
    addParameter(p , 'PermutationTest'        , false, @islogicallike );
    addParameter(p , 'PermutationMethod'      , 'manual', @ischar     );
    addParameter(p , 'PermutationIndex'       , ''   , @ischar        );
    addParameter(p , 'RestrictPermutationByCV', false, @islogicallike );
    % Hyperband (an alternative to grid search for hyperparameter selection)
    addParameter(p , 'HYPERBAND', [] );
    addParameter(p , 'BRACKETS' , [] );
    addParameter(p , 'IterationsPerHyperband', 1000, @isnumeric );
    % Output control
    addParameter(p , 'SmallFootprint', false  , @islogicallike );
    addParameter(p , 'SaveResultsAs'  , 'mat' , @isMatOrJSON   );
    addParameter(p , 'subject_id_fmt' , '%d'  , @ischar        );
    % Debugging (none of which currently work...)
    addParameter(p , 'SanityCheckData'  , []  , @ischar        );
    addParameter(p , 'SanityCheckModel' , []  , @ischar        );
    addParameter(p , 'debug'          , false , @islogicallike );
    % --- searchlight specific --- %
    % N.B. Searchlight is currently broken.
    addParameter(p , 'searchlight'      , 0         , @islogicallike );
    addParameter(p , 'slclassifier'     , ''        , @ischar        );
    addParameter(p , 'slTestToUse'      , ''        , @ischar        );
    addParameter(p , 'slShape'          , ''        , @ischar        );
    addParameter(p , 'slSim_Measure'    , ''        , @ischar        );
    addParameter(p , 'slRadius'         , []        , @isnumeric     );
    addParameter(p , 'slPermutationType', ''        , @ischar        );
    addParameter(p , 'slPermutations'   , 0         , @isscalar      );
    % Parallel only influences the GLMNET operations, and should only be used
    % when running locally. DO NOT USE ON CONDOR.
    addParameter(p , 'PARALLEL'         , false   ,   @islogicallike );
    % Parameters in this section are unused in the analysis, may exist in
    % the parameter file because other progams use them.
    addParameter(p , 'COPY'             , []                         );
    addParameter(p , 'URLS'             , []                         );
    addParameter(p , 'executable'       , []                         );
    addParameter(p , 'wrapper'          , []                         );

    % Parse input parameters ...
    if nargin == 1 && isstruct(varargin{1})
        s = varargin{1};
        x = [fieldnames(s), struct2cell(s)]';
        parse(p, x{:});
    elseif nargin > 0
        % From command line
        parse(p, varargin{:});
    else
        % From json-formatted parameter file
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

    % Check that the correct parameters are passed, given the desired regularization
    regularization = p.Results.regularization;
    switch upper(regularization)
        % The verify_setup_* functions can probably be merged...
        case {'L1L2','GROWL','GROWL2'};
%             assert_required_parameters_RSA(p.Results);
            HYPERPARAMETERS = verify_setup_RSA(regularization, p.Results);
        case {'LASSO','SOSLASSO'};
%             assert_required_parameters_MVPA(p.Results);
            HYPERPARAMETERS = verify_setup_MVPA(regularization, p.Results);
    end

%     DEBUG                   = p.Results.debug;
    PermutationTest         = p.Results.PermutationTest;
    PermutationMethod       = p.Results.PermutationMethod;
    PermutationIndex        = p.Results.PermutationIndex;
%     RestrictPermutationByCV = p.Results.RestrictPermutationByCV;
    SmallFootprint          = p.Results.SmallFootprint;
    RandomSeed              = p.Results.RandomSeed;
    regularization          = p.Results.regularization;
    normalize_wrt            = p.Results.normalize_wrt;
    normalize_data          = p.Results.normalize_data;
    normalize_target        = p.Results.normalize_target;
    BIAS                    = p.Results.bias;
    target_label            = p.Results.target_label;
    target_type             = p.Results.target_type;
    sim_source              = p.Results.sim_source;
    sim_metric              = p.Results.sim_metric;
    filter_labels           = p.Results.filters;
    if iscell(p.Results.data)
        datafiles = p.Results.data;
    else
        datafiles = {p.Results.data};
    end
    data_varname            = p.Results.data_varname;
    cvscheme                = p.Results.cvscheme;
    cvholdout               = p.Results.cvholdout;
    finalholdoutInd         = p.Results.finalholdout;
    orientation             = p.Results.orientation;
    metafile                = p.Results.metadata;
    metadata_varname        = p.Results.metadata_varname;
    tau                     = p.Results.tau;
% --- These are handled in verify_setup* function ---
% --- and stored in HYPERPARAMETERS               ---
%     alpha                   = p.Results.alpha;
%     lambda                  = p.Results.lambda;
%     lamSOS                  = p.Results.alpha;
%     lamL1                   = p.Results.lambda;
%     lambda1                 = p.Results.lambda1;
%     LambdaSeq               = p.Results.LambdaSeq;
    opts                    = p.Results.AdlasOpts;
    SaveResultsAs           = p.Results.SaveResultsAs;
    FMT_subjid              = p.Results.subject_id_fmt;
    % --- searchlight specific --- %
    SEARCHLIGHT        = p.Results.searchlight;
    slSim_Measure      = p.Results.slSim_Measure;
    slPermutationType  = p.Results.slPermutationType;
    slPermutationCount = p.Results.slPermutations;
%     slShape            = p.Results.slShape;
%     slRadius           = p.Results.slRadius;
    % --- HYPERBAND ---
    HYPERBAND = p.Results.HYPERBAND;
    BRACKETS = p.Results.BRACKETS;
    SearchWithHyperband = ~isempty(HYPERBAND);

    if min(RandomSeed) < 1
        RandomSeed = RandomSeed + min(RandomSeed) + 1;
    end
    if ~isempty(RandomSeed) && isscalar(RandomSeed)
        rng(RandomSeed);
    end

    if SEARCHLIGHT && ~strcmpi(slSim_Measure,'nrsa')
        assert(~isempty(slPermutationType));
        assert(~isempty(slPermutationCount));
    end

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
% THIS IS TROUBLING
%    % If cell array with one element, unpack element from cell.
%    datafile = fullfile('D:/MRI/SoundPicture/data/MAT/avg/bysession',uncell(datafile));
%    metafile = fullfile('D:/MRI/SoundPicture/data/MAT/avg/bysession',uncell(metafile));
%
    %% Load metadata
    StagingContainer = load(metafile, metadata_varname);
    metadata = StagingContainer.(metadata_varname); clear StagingContainer;
    subject_label = {metadata.subject};
    [metadata, subjix] = subsetMetadata(metadata, datafiles, FMT_subjid);

    %% Load data   
    X = loadData_new(datafiles, data_varname, FMT_subjid);

    % N.B. Both X and metadata are ordered the same as the 'datafile' cell
    % array. This means the i-th structure in the metadata array corresponds to
    % the i-th cell in the X array. It also means that the order the files were
    % listed dictates the order of these arrays---they are not sorted into
    % ascending numeric or alphabetic order.


    %% Compose and apply filters for each of N subjects
    %%  --- and ---
    %% Load CV indexes, identifying the final holdout set.
    % N.B. the final holdout set is excluded from the rowfilter.
    rowfilter = cell(max(subjix),1);
    colfilter = cell(max(subjix),1);
    cvind     = cell(max(subjix),1);
    cvindAll  = cell(max(subjix),1);
    for i = subjix
        M = selectbyfield(metadata,'subject', subject_label{i});
        if isempty(filter_labels)
            rowfilter{i} = true(1,M.nrow);
            colfilter{i} = true(1,M.ncol);
        else
            [rowfilter{i},colfilter{i}] = composeFilters_new(M.filters, filter_labels);
        end
        cvindAll{i} = M.cvind(:,cvscheme);
        finalholdout = cvindAll{i} == finalholdoutInd;
        % Add the final holdout set to the rowfilter
        rowfilter{i} = forceRowVec(rowfilter{i}) & forceRowVec(~finalholdout);
        % Remove the final holdout set from the cvind, to match.
        cvind{i} = cvindAll{i}(rowfilter{i});
        % Apply the row and column filters
        X{i} = X{i}(rowfilter{i},colfilter{i});
    end

    %% Select targets
    fprintf('\n');
    fprintf('Loading similarity structure\n');
    fprintf('----------------------------\n');
    fprintf('%12s: %s\n', 'target_label', target_label);
    fprintf('%12s: %s\n', 'type', target_type);
    fprintf('%12s: %s\n', 'sim_source', sim_source);
    fprintf('%12s: %s\n', 'sim_metric', sim_metric);
    fprintf('\n');

    tmpS = selectTargets(metadata, target_type, target_label, sim_source, sim_metric, rowfilter(subjix));
    S = cell(max(subjix),1);
    for i = 1:numel(tmpS);
        S(subjix(i)) = tmpS(i);
    end
    clear tmpS;

    % Apply the column filter to the coordinates in the metadata structure
    for i = 1:numel(metadata)
        COORDS = selectbyfield(metadata(i).coords, 'orientation', orientation);
        COORDS_FIELDS = fieldnames(COORDS);
        for j = 1:numel(COORDS_FIELDS)
            cfield = COORDS_FIELDS{j};
            if any(strcmp(cfield, {'ijk','xyz'})) && ~isempty(COORDS.(cfield))
                COORDS.(cfield) = COORDS.(cfield)(colfilter{i},:);
            elseif any(strcmp(cfield, {'ind'})) && ~isempty(COORDS.(cfield))
                COORDS.(cfield) = COORDS.(cfield)(colfilter{i});
            end
        end
        metadata(i).coords = COORDS;
    end
    fprintf('Data Dimensions\n');
    fprintf('%16s%16s%16s\n','subject','initial','filtered');
    fprintf('%s\n',repmat('-',1,16*3));
    for i = subjix
        fprintf('%16d (%6d,%6d) (%6d,%6d)\n',i,numel(rowfilter{i}),numel(colfilter{i}),size(X{i},1),size(X{i},2));
    end
    fprintf('\n');

    fprintf('Data loaded and processed.\n');

    C = cell(max(subjix),1);
    for i = subjix
        switch target_type
            case 'similarity'
                [C{i}, r] = sqrt_truncate_r(S{i}, tau);
                fprintf('S decomposed into %d dimensions (tau=%.2f)\n', r, tau)
            case {'embedding','category'}
                C{i} = S{i};
%                 r = size(C{i},2);
        end
        if size(C{i},1) < size(X{i},1)
            remainder = rem(size(X{i},1), size(C{i},1));
            if remainder > 0
                error('C has fewer rows than X, and number in X is not evenly divisible by number in C.');
            else
                repeatntimes = size(X{i},1) ./ size(C{i},1);
                warning('C has fewer rows than X, and number in X is evenly divisible by number in C. Repeating C %d times to match.', repeatntimes);
                C{i} = repmat(C{i}, repeatntimes, 1);
            end
        end
    end

    % Note on randomization for permutation
    % -------------------------------------
    % A required argument when specifying permutations is a list of
    % "RandomSeeds". These are applied near the beginning of the
    % program (within WholeBrain_RSA), to seed the random number
    % generator.
    %
    % If the PermutationMethod is 'manual', then the RandomSeed has a
    % different (additional) function. It will be used to index into
    % the columns of a n x p matrix, generated in advance, that
    % contains the indexes to generate p unique permutations.
    %
    % In this case, the matrix should stored in a variable named
    % PERMUTATION_INDEXES, contained within a file named
    % PERMUTATION_INDEXES.mat
    fprintf('PermutationTest: %d\n', PermutationTest);
    if PermutationTest
        [a,b] = ndgrid(subjix,RandomSeed);
        Permutations = struct('subject', num2cell(a),'RandomSeed', num2cell(b),'index',[]);
        switch PermutationMethod
%                 case 'simple'
%                     if RestrictPermutationByCV
%                         C = permute_target(C, PermutationMethod, cvind);
%                     else
%                         C = permute_target(C, PermutationMethod);
%                     end
            case 'manual'
                StagingArea = load(PermutationIndex, 'PERMUTATION_INDEX');
                PERMUTATION_INDEX = StagingArea.PERMUTATION_INDEX;
                for i = subjix
                    for j = RandomSeed
                        P = selectbyfield(Permutations,'subject',i,'RandomSeed',j);
                        %  This is kind of a hack to handle the fact eliminating
                        %  outlying rows and rows belonging to the final holdout
                        %  set will create gaps in the index.
                        [~, ix] = sort(PERMUTATION_INDEX{i}(rowfilter{i}, j));
                        [~, permutation_index] = sort(ix);
                        if size(permutation_index, 1) < size(C{i}, 1)
                            permutation_index = extend_permutation_index(permutation_index, size(C{i}, 1));
                        end
                        P.index = permutation_index;
                        Permutations = replacebyfield(Permutations, P, 'subject', i, 'RandomSeed', j);
                    end
                end
            otherwise
                error('crcox:NotImplemented', 'Permutations need to be specified manually.');
        end
    else
        RandomSeed = 0;
        Permutations = struct('subject', num2cell(subjix), 'RandomSeed', 0, 'index',[]);
        for i = subjix
            P = selectbyfield(Permutations,'subject',i,'RandomSeed',0);
            P.index = (1:size(C{i}, 1))';
            Permutations = replacebyfield(Permutations, P, 'subject', i, 'RandomSeed', 0);
        end
    end

    %% --- Setting regularization parameters and running models ---
    % This is being handled within the 
    switch upper(regularization)
        case {'GROWL','GROWL2','L1L2','LASSO'}
            SubjectsParameter = subjix;
        case 'SOSLASSO'
            SubjectsParameter = {{subjix}};
    end

    if exist('checkpoint.mat','file')
        load('checkpoint.mat', 'ModelInstances', 'bracket_index');
    else
        ModelInstances = ModelContainer( ...
            'subject'          , SubjectsParameter, ...
            'RandomSeed'       , RandomSeed       , ...
            'cvholdout'        , cvholdout        , ...
            'finalholdout'     , finalholdoutInd  , ...
            'bias'             , BIAS             , ...
            'target_label'     , target_label     , ...
            'target_type'      , target_type      , ...
            'sim_metric'       , sim_metric       , ...
            'sim_source'       , sim_source       , ...
            'normalize_data'   , normalize_data   , ...
            'normalize_target' , normalize_target , ...
            'normalize_wrt'    , normalize_wrt    , ...
            'regularization'   , regularization   , ...
            'HYPERPARAMETERS'  , HYPERPARAMETERS);
        % TODO: There is probably a smart way to incorporate this
        % functionality (basically, child fields that are associated with a
        % parent field) within the ModelContainer expansion function.
        for i = 1:numel(ModelInstances)
            if iscell(datafiles) && numel(datafiles) > 1;
                ModelInstances(i).data = datafiles(ModelInstances(i).subject);
            else
                ModelInstances(i).data = ascell(datafiles);
            end
            if iscell(data_varname) && numel(data_varname) > 1;
                ModelInstances(i).data_varname = data_varname(ModelInstances(i).subject);
            else
                ModelInstances(i).data_varname = ascell(data_varname);
            end
            if iscell(metafile) && numel(metafile) > 1;
                ModelInstances(i).metadata = metafile(ModelInstances(i).subject);
            else
                ModelInstances(i).metadata = ascell(metafile);
            end
            if iscell(metadata_varname) && numel(metadata_varname) > 1;
                ModelInstances(i).metadata_varname = metadata_varname(ModelInstances(i).subject);
            else
                ModelInstances(i).metadata_varname = ascell(metadata_varname);
            end
        end
        bracket_index = 1;
    end

    xyz = cell(numel(metadata),1);
    for ii = 1:numel(metadata)
        z = strcmp({metadata(ii).coords.orientation}, orientation);
        xyz{ii} = metadata(ii).coords(z).xyz;
    end
    switch upper(regularization)
        case 'SOSLASSO'
            for ii = 1:numel(ModelInstances)
                diameter = ModelInstances(ii).diameter;
                overlap = ModelInstances(ii).overlap;
                shape = ModelInstances(ii).shape;
                ModelInstances(ii).G = coordGrouping(xyz, diameter, overlap, shape);
            end

        case 'LASSO'
            for ii = 1:numel(ModelInstances)
                ModelInstances(ii).G = coordGrouping(xyz, 0, 0, 'unitary');
            end
    end

    if SearchWithHyperband
        n = BRACKETS.n;
        r = BRACKETS.r;
        while 1
            opts.max_iter = r(bracket_index) * p.Results.IterationsPerHyperband; % This 1000 is an important constant... might want to think about this/expose it as a parameter.
            ModelInstances = learn_similarity_encoding(ModelInstances, C, X, regularization,...
                'cvind'          , cvind        , ...
                'permutations'   , Permutations , ...
                'AdlasOpts'      , opts);
            % Delete low ranked configurations:
            ModelInstances = hyperband_pick_top_n(ModelInstances, n(bracket_index));
            bracket_index = bracket_index + 1;
            if bracket_index < numel(n)
                save('checkpoint.mat', 'ModelInstances', 'bracket_index');
            else
                if exist('checkpoint.mat', 'file')
                    delete('checkpoint.mat');
                end
                break
            end
        end
    else
        % Grid search
        ModelInstances = learn_similarity_encoding(ModelInstances, C, X, regularization,...
            'cvind'          , cvind        , ...
            'permutations'   , Permutations , ...
            'AdlasOpts'      , opts);
    end

    cur = 0;
    n = numel(ModelInstances);
    for i = 1:n
        MODEL = ModelInstances(i).Model;
        CONTEXT = rmfield(ModelInstances(i), 'Model');
        if i == 1
            % Preallocate on first pass
            [results,nResultsPerModel] = MODEL.getResults(CONTEXT,metadata,'Initialize',n);
        end
        a = cur + 1;
        b = cur + nResultsPerModel;
        results(a:b) = MODEL.getResults(CONTEXT,metadata);
        cur = a;
    end
    %% --- Package results ---
%     results = repmat(struct( ...
%         'Uz'               , [] , ...
%         'Cz'               , [] , ...
%         'Sz'               , [] , ...
%         'target_label'     , [] , ...
%         'target_type'      , [] , ...
%         'sim_source'       , [] , ...
%         'sim_metric'       , [] , ...
%         'data'             , [] , ...
%         'data_varname'     , [] , ...
%         'metadata'         , [] , ...
%         'metadata_varname' , [] , ...
%         'subject'          , [] , ...
%         'cvholdout'        , [] , ...
%         'finalholdout'     , [] , ...
%         'regularization'   , [] , ...
%         'lambda'           , [] , ...
%         'lambda1'          , [] , ...
%         'LambdaSeq'        , [] , ...
%         'tau'              , [] , ...
%         'bias'             , [] , ...
%         'normalize_wrt'     , [] , ...
%         'normalize_data'   , [] , ...
%         'normalize_target' , [] , ...
%         'nz_rows'          , [] , ...
%         'nzv'              , [] , ...
%         'nvox'             , [] , ...
%         'coords'           , [] , ...
%         'err1'             , [] , ...
%         'err2'             , [] , ...
%         'iter'             , [] ), numel(ModelInstances), 1);
%     for iResult = 1:numel(ModelInstances)
%         A = ModelInstances(iResult);
%         if A.bias
%             Uz = A.Model.X(1:end-1,:);
%         else
%             Uz = A.Model.X;
%         end
%         if ~SmallFootprint
%             results(iResult).coords = COORDS;
%             ix = find(any(Uz, 2));
%             for j = 1:numel(COORDS_FIELDS)
%                 cfield = COORDS_FIELDS{j};
%                 if any(strcmp(cfield, {'ijk','xyz'})) && ~isempty(COORDS.(cfield))
%                     results(iResult).coords.(cfield) = COORDS.(cfield)(ix,:);
%                 elseif any(strcmp(cfield, {'ind'})) && ~isempty(COORDS.(cfield))
%                     results(iResult).coords.(cfield) = COORDS.(cfield)(ix);
%                 end
%             end
%             results(iResult).Uz = Uz;
%             results(iResult).Cz = A.Model.A * A.Model.X;
%         end
%         results(iResult).subject = A.subject;
%         results(iResult).bias = A.bias;
%         results(iResult).nz_rows = any(Uz,2);
%         results(iResult).nzv = nnz(results(iResult).nz_rows);
%         results(iResult).nvox = numel(results(iResult).nz_rows);
%         results(iResult).cvholdout = A.cvholdout;
%         results(iResult).finalholdout = finalholdoutInd;
%         results(iResult).lambda = A.lambda;
%         results(iResult).lambda1 = A.lambda1;
%         results(iResult).LambdaSeq = A.LambdaSeq;
%         results(iResult).regularization = A.regularization;
%         results(iResult).tau = tau;
%         results(iResult).normalize_data = A.normalize_data;
%         results(iResult).normalize_target = normalize_target;
%         results(iResult).normalize_wrt = A.normalize_wrt;
%         results(iResult).data = p.Results.data;
%         results(iResult).data_var = p.Results.data_varname;
%         results(iResult).metadata = p.Results.metadata;
%         results(iResult).metadata_var = p.Results.metadata_varname;
%         results(iResult).target_label = p.Results.target_label;
%         results(iResult).target_type = p.Results.target_type;
%         results(iResult).sim_source = p.Results.sim_source;
%         results(iResult).sim_metric = p.Results.sim_metric;
%         results(iResult).err1 = A.Model.testError;
%         results(iResult).err2 = A.Model.trainingError;
%         results(iResult).iter = A.Model.iter;
%         results(iResult).RandomSeed = A.RandomSeed;
% %         results(iResult).RandomSeed = p.Results.PermutationIndex;
%     end

    fprintf('Saving stuff...\n');

    %% Save results
    rinfo = whos('results');
    switch SaveResultsAs
        case 'mat'
            if rinfo.bytes > 2e+9 % 2 GB
                save('results.mat','results','-v7.3');
            else
                save('results.mat','results');
            end
        case 'json'
            if rinfo.bytes > 16e+6 % 16 MB
                disp('WARNING: Results structure too large to save as JSON (excedes MongoDB 16MB limit). Saving as .mat...')
                if rinfo.bytes > 2e+9 % 2 GB
                    save('results.mat','results','-v7.3');
                else
                    save('results.mat','results');
                end
            else
                savejson('',results,'FileName','results.json','ForceRootName',false);
            end
    end

    fprintf('Done!\n');

%     function condition_handling_searchlight()
%         X = uncell(X);
%         S = uncell(S);
%         cvind = uncell(cvind);
%         cvset = unique(cvind);
%         colfilter = uncell(colfilter);
% 
%         % create a 3D binary mask
%         [mask,dxyz] = coordsTo3dMask(metadata.coords.xyz);
% 
%         % Translate slradius (in mm) to sl voxels
%         % N.B. Because voxels need not be symmetric cubes, but Seachmight will
%         % generate symmetric spheres from a single radius parameter, we need to
%         % select one value of the three that will be produced in this step. I am
%         % arbitrarily choosing the max, to err on the side of being inclusive.
%         slradius_ijk = max(round(slRadius ./ dxyz));
% 
%         % create the "meta" neighbourhood structure
%         meta = createMetaFromMask(mask, 'radius', slradius_ijk);
%         labels = metadata.itemindex(rowfilter);
%         labelsRun = metadata.runindex(rowfilter);
% 
%         results.similarity_measure = slSim_Measure;
%         if strcmpi('nrsa',slSim_Measure)
%             error('Searchlight Network RSA is not implemented properly yet. Exiting...');
% 
%         else
%             fprintf('PermutationTest: %d\n', PermutationTest);
%             if PermutationTest
%                 for ic = unique(cvind)'
%                     fprintf('Permuting CV %d...\n', ic);
%                     s = S(cvind==ic, cvind==ic);
%                     n = size(s,1);
%                     permix = randperm(n);
%                     S(cvind==ic, cvind==ic) = S(permix, permix);
%                 end
%             end
% 
%             [structureScoreMap] = computeSimilarityStructureMap(...
%                 slSim_Measure,...
%                 X,labels,...
%                 X,labels,...
%                 'meta',meta,'similarityStructure',S,...
%                 'permutationTest',slPermutationType, slPermutationCount,...
%                 'groupLabels',labelsRun,labelsRun);
% 
%             results.structureScoreMap = structureScoreMap;
%             results.RandomSeed = RandomSeed;
%         end
% 
%         for i_nested = 1:numel(results)
%             results(i_nested).coords = COORDS;
%         end
%     end
end

function [hyperparameters] = verify_setup_MVPA(regularization, p)
    % Each regularization requires different lambda configurations. This private
    % function ensures that everything has been properly specified.
    SearchWithHyperband = ~isempty(p.HYPERBAND);
    switch upper(regularization)
        case 'LASSO_GLMNET'
            if isfield(p,'lambda') && isempty(p.lambda)
                warning('Lamba was not specified. GLMnet will attempt to determine lambda1 through cross validation.');
                lambda = nan(1);
            else
                lambda = p.lambda;
            end
            if isfield(p,'alpha') && ~isempty(p.alpha) && (p.alpha ~= 1)
                warning('GLMnet performs lasso when alpha=1. Forcing alpha=1.');
                alpha = 1;
            else
                alpha = p.alpha;
            end
            hyperparameters = struct('alpha',alpha,'lambda',lambda,'hyperband',SearchWithHyperband);

        case 'LASSO'
            assert(any(isfield(p,{'lambda','lamL1'})), 'Either lambda or lamL1 (synonymns for LASSO regularization) must be defined.');
            if isfield(p,'alpha') && ~isempty(p.alpha)
                warning('Alpha is not relevant for performing Lasso with SOSLasso_logistic. Forcing lamSOS=1, which means that lamL1 will not be scaled up or down. The critical thing for lasso is that all voxels get their own group. This is enforced elsewhere...');
                lamSOS = 0;
                lamL1 = p.lambda;
            end
            if isfield(p,'lamL1') && ~isempty(p.lamL1)
                if isfield(p,'lambda') && ~isempty(p.lambda)
                    warning('lamL1 takes precedence over lambda when both are defined. Lambda is being ignored.');
                end
                lamL1 = p.lamL1;
            elseif isfield(p,'lambda') && ~isempty(p.lambda)
                lamL1 = p.lambda;
            end
            if isfield(p,'lamSOS') && ~isempty(p.lamSOS)
                if p.lamSOS ~= 1
                    warning('lamSOS is not, strictly speaking, relevant for performing Lasso with SOSLasso_logistic. It will act as a scalar modifier on lamL1, which because lamL1 is just a constant anyway, is pointless. Forcing lamSOS=1 to prevent scaling lamL1.');
                    lamSOS = 1;
                else
                    lamSOS = p.lamSOS;
                end
            end
            hyperparameters = struct('lamSOS',lamSOS,'lamL1',lamL1,'hyperband',SearchWithHyperband);

        case 'SOSLASSO'
            assert(all(isfield(p,{'diameter','shape','overlap'})) && ~isempty(p.diameter) && ~isempty(p.shape) && ~isempty(p.overlap), 'diameter, shape, and overlap must all be defined for SOSLASSO regularization.');
            assert( ...
                all(isfield(p,{'alpha','lambda'})) && ~isempty(p.alpha) && ~isempty(p.lambda) || ...
                all(isfield(p,{'lamSOS','lamL1'})) && ~isempty(p.lamSOS) && ~isempty(p.lamL1), ...
                'Either (alpha and lambda) or (lamSOS and lamL1) must be defined for SOSLASSO regularization.');
            if any(isfield(p,{'alpha','lambda'})) && ~all(isfield(p,{'lamSOS','lamL1'}))
                assert(all(isfield(p,{'alpha','lambda'})), 'Alpha and lambda must both be specified to generate lamSOS and lamL1 for SOSLASSO regularization.')
            end
            if all(isfield(p,{'alpha','lambda'})) && ~isempty(p.alpha) && ~isempty(p.lambda)
                warning('Alpha will be translated to lamSOS with respect to lambda.');
                warning('Lambda will be translated to lamL1 with respect to alpha.');
                [lamSOS,lamL1] = ratio2independent(p.alpha,p.lambda);
            elseif all(isfield(p,{'lamSOS','lamL1'})) && ~isempty(p.lamSOS) && ~isempty(p.lamL1)
                lamSOS = p.lamSOS;
                lamL1 = p.lamL1;
            end
            hyperparameters = struct('lamSOS',lamSOS,'lamL1',lamL1,'diameter',p.diameter,'shape',p.shape,'overlap',p.overlap,'hyperband',SearchWithHyperband);

    end
end

function [hyperparameters] = verify_setup_RSA(regularization, p)
    % Each regularization requires different lambda configurations. This private
    % function ensures that everything has been properly specified.
    SearchWithHyperband = ~isempty(p.HYPERBAND);
    switch upper(regularization)
        case 'NONE'
            if ~isempty(lambda) || ~isempty(lambda1)
                warning('Regularization was set to none, but lambda values were provided. They will be ignored.')
            end
            hyperparameters = struct('hyperband',false);

        case 'L1L2_GLMNET'
            if isfield(p,'lambda') && ~isempty(p.lambda)
                lam = p.lambda;
            else
                warning('Lamba was not specified. GLMnet will attempt to determine lambda1 through cross validation.');
                lam = nan(1);
            end
            hyperparameters = struct('alpha',1,'lambda',lam,'hyperband',SearchWithHyperband);

        case 'L1L2'
            if isfield(p,'lambda1') && ~isempty(p.lambda1)
                warning('Group Lasso does not use the lambda1 parameter. It is being ignored.');
            end
            if isfield(p,'LambdaSeq') && ~isempty(p.LambdaSeq) && all(~stcmp(p.LambdaSeq,'none'))
                warning('Group Lasso does not use a lambda sequence. Setting to ''none''.');
            end
            assert(~isempty(p.lambda), 'Group Lasso requires lambda.');
            lam    = p.lambda;
            lam1   = [];
            lamSeq = 'none';
            hyperparameters = struct('lambda',lam,'lambda1',lam1,'lambdaSeq',lamSeq,'hyperband',SearchWithHyperband);

        case {'GROWL','GROWL2'}
            assert(isfield(p,'lambda') && ~isempty(p.lambda) && any(~isnan(p.lambda)), 'grOWL requires lambda.');
            assert(isfield(p,'lambda1') && ~isempty(p.lambda1) && any(~isnan(p.lambda1)), 'grOWL requires lambda1.');
            assert(isfield(p,'LambdaSeq') && ~isempty(p.LambdaSeq), 'A LambdaSeq type (linear or exponential) must be set when using grOWL*.');
            lam    = p.lambda;
            lam1   = p.lambda1;
            lamSeq = p.LambdaSeq;
            hyperparameters = struct('lambda',lam,'lambda1',lam1,'lambdaSeq',lamSeq,'hyperband',SearchWithHyperband);

    end
end

% function assertRequiredParameters_RSA(params)
%     required = {'target_label','sim_metric','sim_source','data', ...
%         'metadata','cvscheme','cvholdout','finalholdout','orientation'};
%     N = length(required);
%     for i = 1:N
%         req = required{i};
%         assert(isfield(params,req), '%s must exist in params structure! Exiting.',req);
%         assert(~isempty(params.(req)), '%s must be set. Exiting.',req);
%     end
% end

function permutation_index = extend_permutation_index(permutation_index, target_length)
    remainder = rem(target_length, size(permutation_index, 1));
    if remainder > 0
        error('permutation_index has fewer rows than C, and number in C is not evenly divisible by number in permutation_index.');
    else
        repeatntimes = target_length / size(permutation_index,1);
        warning('permutation_index has fewer rows than C, and number in C is evenly divisible by number in permutation_index. Repeating permutation_index %d times to match.', repeatntimes);
        CC = cell(repeatntimes, 1);
        for k = 1:repeatntimes
            CC{k} = permutation_index + (size(permutation_index, 1) * (k-1));
        end
        permutation_index = cell2mat(CC);
    end
end
function b = islogicallike(x)
    b = any(x == [1,0]);
end

function b = isintegerlike(x)
    b = mod(x,1) == 0;
end

function b = isMatOrJSON(x)
    b = any(strcmpi(x, {'mat','json'}));
end
% function b = isscalarOrEmpty(x)
%     b = isscalar(x) || isempty(x);
% end

% OTHER REGULARIZATION METHODS
% ============================
% case 'L1L2_GLMNET'
%     if isempty(gcp('nocreate')) && PARALLEL
%         ppp = parpool('local');
%     end
%     [results,info] = learn_similarity_encoding(S, X, regularization, target_type,...
%         'tau'            , tau            , ...
%         'lambda1'        , lambda1        , ...
%         'cvind'          , cvind          , ...
%         'cvholdout'      , cvholdout      , ...
%         'normalize'      , normalize      , ...
%         'bias'           , BIAS           , ...
%         'DEBUG'          , DEBUG          , ...
%         'SmallFootprint' , SmallFootprint , ...
%         'AdlasOpts'      , opts); %#ok<ASGLU>
%     if ~isempty(gcp('nocreate')) && PARALLEL && (exist('ppp', 'var') == 1)
%         delete(ppp);
%     end

% OLD RESULTS SETUP AND ANALYSIS LAUNCHING
% ========================================
%             % Define results structure
%             results.Uz = [];
%             results.Cz = [];
%             results.Sz = [];
%             results.nz_rows =  [];
%             results.target_label = target_label;
%             results.subject =  [];
%             results.cvholdout = [];
%             results.finalholdout = [];
%             results.lambda = [];
%             results.lambda1 = [];
%             results.LambdaSeq = [];
%             results.regularization = [];
%             results.bias = [];
%             results.normalize = [];
%             results.nzv = [];
%             %      results.p1      =  [];
%             %      results.p2      =  [];
%             %      results.cor1    =  [];
%             %      results.cor2    =  [];
%             %      results.p1t     =  [];
%             %      results.p2t     =  [];
%             %      results.cor1t   =  [];
%             %      results.cor2t   =  [];
%             results.coords  = [];
%             results.structureScoreMap = zeros(1, size(meta.voxelsToNeighbours,1));
%             results.structurePvalueMap = zeros(1, size(meta.voxelsToNeighbours,1));
%             results.err1    =  zeros(1, size(meta.voxelsToNeighbours,1));
%             results.err2    =  zeros(1, size(meta.voxelsToNeighbours,1));
%             results.iter    =  [];
% 
%             % Preallocate
%             if isempty(lambda); nlam = 1; else nlam = numel(lamba); end
%             if isempty(lambda1); nlam1 = 1; else nlam1 = numel(lambda1); end
%             results(numel(cvset)*nlam*nlam1).Uz = [];
% 
%             for iVolume = 1:size(meta.voxelsToNeighbours,1)
%                 sl = meta.voxelsToNeighbours(iVolume,1:meta.numberOfNeighbours(iVolume));
%                 switch upper(regularization)
% %                     case 'L1L2'
% %                         [lambda1, err_L1L2] = fminbnd(@(x) optimizeGroupLasso(S,X(:,sl),tau,cvind,cvholdout,normalize,PermutationTest,x), 0, 32);
% %                         if finalholdout > 0
% %                             [tmpr,info] = learn_similarity_encoding(Sall, Xall(:,sl), regularization, target_type,...
% %                                 'tau'            , tau            , ...
% %                                 'lambda1'        , lambda1        , ...
% %                                 'cvind'          , cvindAll       , ...
% %                                 'cvholdout'      , finalholdoutInd, ...
% %                                 'normalize'      , normalize      , ...
% %                                 'bias'           , BIAS           , ...
% %                                 'DEBUG'          , DEBUG          , ...
% %                                 'PermutationTest', PermutationTest, ...
% %                                 'PermutationMethod', PermutationMethod, ...
% %                                 'RestrictPermutationByCV', RestrictPermutationByCV, ...
% %                                 'SmallFootprint' , SmallFootprint , ...
% %                                 'AdlasOpts'      , opts); %#ok<ASGLU>
% %                         else
% %                             [tmpr,info] = learn_similarity_encoding(S, X(:,sl), regularization, target_type,...
% %                                 'tau'            , tau            , ...
% %                                 'lambda1'        , lambda1        , ...
% %                                 'cvind'          , cvind          , ...
% %                                 'cvholdout'      , cvholdout      , ...
% %                                 'normalize'      , normalize      , ...
% %                                 'bias'           , BIAS           , ...
% %                                 'DEBUG'          , DEBUG          , ...
% %                                 'PermutationTest', PermutationTest, ...
% %                                 'PermutationMethod', PermutationMethod, ...
% %                                 'RestrictPermutationByCV', RestrictPermutationByCV, ...
% %                                 'SmallFootprint' , SmallFootprint , ...
% %                                 'AdlasOpts'      , opts); %#ok<ASGLU>
% %                         end
% 
%                     case 'GROWL'
%                         [tmpr,info] = learn_similarity_encoding(C, X(:,sl), regularization, target_type,...
%                             'tau'            , tau            , ...
%                             'lambda'         , lambda         , ...
%                             'LambdaSeq'      , LambdaSeq      , ...
%                             'cvind'          , cvind          , ...
%                             'cvholdout'      , cvholdout      , ...
%                             'normalize'      , normalize      , ...
%                             'bias'           , BIAS           , ...
%                             'DEBUG'          , DEBUG          , ...
%                             'SmallFootprint' , SmallFootprint , ...
%                             'AdlasOpts'      , opts); %#ok<ASGLU>
% 
%                     case 'GROWL2'
%                         [tmpr,info] = learn_similarity_encoding(C, X(:,sl), regularization, target_type,...
%                             'tau'            , tau            , ...
%                             'lambda'         , lambda         , ...
%                             'lambda1'        , lambda1        , ...
%                             'LambdaSeq'      , LambdaSeq      , ...
%                             'cvind'          , cvind          , ...
%                             'cvholdout'      , cvholdout      , ...
%                             'normalize'      , normalize      , ...
%                             'bias'           , BIAS           , ...
%                             'DEBUG'          , DEBUG          , ...
%                             'SmallFootprint' , SmallFootprint , ...
%                             'AdlasOpts'      , opts); %#ok<ASGLU>
%                 end
%                 for iResult = 1:numel(tmpr)
%                     results(iResult).err1(iVolume) = tmpr(iResult).err1;
%                     results(iResult).err2(iVolume) = tmpr(iResult).err2;
%                     results(iResult).structureScoreMap(iVolume) = tmpr(iResult).structureScoreMap;
%                 end
%             end
