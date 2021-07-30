function WISC_MVPA(varargin)
    p = inputParser;
    p.KeepUnmatched = false;
    % ----------------------Set parameters---------------------------------
    % Model definition
    addParameter(p , 'regularization'   , []      , @ischar        );
    addParameter(p , 'bias'             , false   , @islogicallike );
    addParameter(p , 'alpha'            , []      , @isnumeric     );
    addParameter(p , 'lambda'           , []      , @isnumeric     );
    addParameter(p , 'lambda1'          , []      , @isnumeric     );
    addParameter(p , 'lamSOS'           , []      , @isnumeric     );
    addParameter(p , 'lamL1'            , []      , @isnumeric     );
    addParameter(p , 'lamL2'            , []      , @isnumeric     );
    addParameter(p , 'LambdaSeq'        , []      , @ischar        );
    addParameter(p , 'AdlasOpts'        , struct(), @isstruct      );
    addParameter(p , 'diameter'         , []      , @isnumeric     );
    addParameter(p , 'shape'            , []      , @ischar        );
    addParameter(p , 'overlap'          , []      , @isnumeric     );
    addParameter(p , 'sosgroups'        , []      , @iscellstr     );
    % Target definition
    addParameter(p , 'target_label'     , [], @ischar );
    addParameter(p , 'target_type'      , [], @ischar );
    addParameter(p , 'sim_source'       , [], @ischar );
    addParameter(p , 'sim_metric'       , [], @ischar );
    addParameter(p , 'tau'              , [], @isnumeric     );
    addParameter(p , 'FiltersToApplyBeforeEmbedding' , [], @(x) ischar(x) || iscellstr(x) || isstring(x));
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
    addParameter(p , 'normalize_data'         , 'none'         , @ischar );
    addParameter(p , 'normalize_target'       , 'none'         , @ischar );
    addParameter(p , 'normalize_wrt'          , 'all_examples' , @ischar );
    addParameter(p , 'scale_singular_vectors' , true , @islogicallike ); % NEW VARIABLE!
    % Permutation
    addParameter(p , 'RandomSeed'             , 0, @(x) isnumeric(x) && all(x>=0));
    addParameter(p , 'PermutationTest'        , false, @islogicallike );
    addParameter(p , 'PermutationMethod'      , 'manual' , @ischar    );
    addParameter(p , 'PermutationIndex'       , ''       , @ischar    );
    addParameter(p , 'perm_varname' , 'PERMUTATION_INDEX', @ischar    );
    addParameter(p , 'RestrictPermutationByCV', false, @islogicallike );
    % Hyperband (an alternative to grid search for hyperparameter selection)
    addParameter(p , 'SearchWithHyperband', false );
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
    addParameter(p , 'HYPERBAND'  , [] );
    addParameter(p , 'COPY'       , [] );
    addParameter(p , 'URLS'       , [] );
    addParameter(p , 'executable' , [] );
    addParameter(p , 'wrapper'    , [] );

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

    % --- setup hyperparameters ---
    switch upper(p.Results.regularization)
        % The verify_setup_* functions can probably be merged...
        case {'L1L2','GROWL','GROWL2'}
            HYPERPARAMETERS = verify_setup_RSA(p.Results);
        case {'RIDGE','LASSO','SOSLASSO'}
            HYPERPARAMETERS = verify_setup_MVPA(p.Results);
    end
    % --- searchlight specific ---
    if p.Results.searchlight && ~strcmpi(p.Results.slSim_Measure,'nrsa')
        assert(~isempty(p.Results.slPermutationType));
        assert(~isempty(p.Results.slPermutations));
    end
    % --- Initialize the random number generator, as needed ---
    if ~isempty(p.Results.RandomSeed) && isscalar(p.Results.RandomSeed) && ~strcmpi(p.Results.PermutationMethod,'manual')
        rng(p.Results.RandomSeed);
    end
    % If values originated in a YAML file, and scientific notation is used,
    % the value may have been parsed as a string. Check and correct.
    opts = p.Results.AdlasOpts;
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

    %% Load data
    % Load metadata object 
    StagingContainer = load(p.Results.metadata, p.Results.metadata_varname);
    metadata = StagingContainer.(p.Results.metadata_varname);
    clear StagingContainer;
    % Load permutation object
    if p.Results.PermutationTest && strcmpi(p.Results.PermutationMethod,'manual')
    	StagingContainer = load(p.Results.PermutationIndex, p.Results.perm_varname);
        permutations = StagingContainer.(p.Results.perm_varname);
        clear StagingContainer;
    end
    % Loop over datafiles
    datafiles  = ascell(p.Results.data);
    N = length(datafiles);
    SubjectArray = repmat(Subject, 1, N);
    for i = 1:N
        x = datafiles{i};
        fprintf('Loading %s from  %s...\n', p.Results.data_varname, x);
        % Log and process data filename and data label
        SubjectArray(i) = SubjectArray(i).setFilename(datafiles{i});
        SubjectArray(i) = SubjectArray(i).setLabel(p.Results.data_varname);
        SubjectArray(i) = SubjectArray(i).setDataFromFilenameAndLabel();
        SubjectArray(i) = SubjectArray(i).setSubjectIDFromFilename(p.Results.subject_id_fmt);
        % Select metadata matching current subject ID
        M = selectbyfield(metadata, 'subject', SubjectArray(i).subject);
        if numel(M) > 1
            error('WISC_MVPA:SubjectSelect','Subject ID %s does not uniquely match a metadata entry.', num2str(SubjectArray(i).subject));
        end
        % Pull content from metadata object
        CV = M.cvind(:,p.Results.cvscheme);
        RF = selectbyfield(M.filters, 'label', [p.Results.filters,p.Results.FiltersToApplyBeforeEmbedding], 'dimension', 1);
        CF = selectbyfield(M.filters, 'label', p.Results.filters, 'dimension', 2);
        T = selectbyfield(M.targets, ...
            'label', p.Results.target_label, ...
            'type', p.Results.target_type, ...
            'sim_source', p.Results.sim_source, ...
            'sim_metric', p.Results.sim_metric);
        % Set metadata for subject
        SubjectArray(i) = SubjectArray(i).setCVScheme(CV);
        SubjectArray(i) = SubjectArray(i).setRowFilters(RF);
        SubjectArray(i) = SubjectArray(i).setColFilters(CF);
        SubjectArray(i) = SubjectArray(i).setFinalHoldoutFilter(p.Results.finalholdout);
        SubjectArray(i) = SubjectArray(i).setCoords(M.coords);
        SubjectArray(i) = SubjectArray(i).setTargets(T);
        disp(p.Results.scale_singular_vectors)
        if ~isempty(p.Results.tau) && p.Results.tau ~= 0
            % If tau is set, generate embeddings from target similarity matrices.
            if isempty(p.Results.FiltersToApplyBeforeEmbedding)
                SubjectArray(i) = SubjectArray(i).generateEmbeddings(p.Results.tau, 'ApplySingularValues', p.Results.scale_singular_vectors);
            else
                SubjectArray(i) = SubjectArray(i).generateEmbeddings(p.Results.tau, 'ApplySingularValues', p.Results.scale_singular_vectors, 'PreFilter', p.Results.FiltersToApplyBeforeEmbedding);
            end
        end
        % Remove Pre-filters if they are not also used as standard filters.
        for j = 1:numel(p.Results.FiltersToApplyBeforeEmbedding)
            xx = p.Results.FiltersToApplyBeforeEmbedding{j};
            if ~strcmp(xx, p.Results.filters)
                z = strcmp(xx,{SubjectArray(i).rowfilters.label});
                SubjectArray(i).rowfilters(z) = [];
            end
        end
        % Set permutation data. If a permutation test is not being run,
        % RandomSeed == 0 and a dummy permutation index will be encoded
        % that is simply 1:Subject.nTotalExamples.
        if p.Results.PermutationTest && strcmpi(p.Results.PermutationMethod, 'manual')
            P = selectbyfield(permutations, 'subject', SubjectArray(i).subject);
            SubjectArray(i) = SubjectArray(i).setPermutations(p.Results.PermutationMethod, p.Results.RandomSeed, P.permutation_index);
        else
            SubjectArray(i) = SubjectArray(i).setPermutations('none', p.Results.RandomSeed);
        end
    end

    % Report target infomation
    LogicalString = {'FALSE', 'TRUE'};
    fprintf('\n');
    fprintf('Target Structure Summary\n');
    fprintf('------------------------\n');
    fprintf('%24s: %s\n', 'target_label', p.Results.target_label);
    fprintf('%24s: %s\n', 'type', p.Results.target_type);
    fprintf('%24s: %s\n', 'sim_source', p.Results.sim_source);
    fprintf('%24s: %s\n', 'sim_metric', p.Results.sim_metric);
    fprintf('%24s: %s\n', 'scale_singular_vectors', LogicalString{p.Results.scale_singular_vectors + 1});
    fprintf('\n');

    % Report Data information
    fprintf('Data Dimensions\n');
    fprintf('---------------\n');
    fprintf('%16s%16s%16s\n','subject','initial','filtered');
    fprintf('%s\n',repmat('-',1,16*3));
    for i = 1:numel(SubjectArray)
        fprintf('%16s (%6d,%6d) (%6d,%6d)\n',num2str(SubjectArray(i).subject),size(SubjectArray(i).getData('unfiltered',true)),size(SubjectArray(i).getData('unfiltered',false)));
    end
    fprintf('\n');
    fprintf('Data loaded and processed.\n');

    %% --- Setting regularization parameters and running models ---
    % This is being handled within the 
    switch upper(p.Results.regularization)
        case {'GROWL','GROWL2','L1L2','LASSO','RIDGE'}
            SubjectsParameter = [SubjectArray.subject];
        case 'SOSLASSO'
            SubjectsParameter = 1;
    end

    if exist('checkpoint.mat','file')
        load('checkpoint.mat', 'ModelInstances', 'bracket_index');
    else
        ModelInstances = ModelContainer( ...
            'subject'                , SubjectsParameter                , ...
            'RandomSeed'             , p.Results.RandomSeed             , ...
            'cvholdout'              , p.Results.cvholdout              , ...
            'finalholdout'           , p.Results.finalholdout           , ...
            'bias'                   , p.Results.bias                   , ...
            'target_label'           , p.Results.target_label           , ...
            'target_type'            , p.Results.target_type            , ...
            'sim_metric'             , p.Results.sim_metric             , ...
            'sim_source'             , p.Results.sim_source             , ...
            'normalize_data'         , p.Results.normalize_data         , ...
            'normalize_target'       , p.Results.normalize_target       , ...
            'normalize_wrt'          , p.Results.normalize_wrt          , ...
            'scale_singular_vectors' , p.Results.scale_singular_vectors , ...
            'regularization'         , p.Results.regularization         , ...
            'HYPERPARAMETERS'        , HYPERPARAMETERS);
        % TODO: There is probably a smart way to incorporate this
        % functionality (basically, child fields that are associated with a
        % parent field) within the ModelContainer expansion function.
        if strcmpi(p.Results.regularization, 'SOSLASSO')
            % The assumption is that, when running SOSLASSO, all subjects
            % are associated with each model instance.
            [ModelInstances.data] = deal({SubjectArray.filename});
        else
            for i = 1:numel(SubjectArray)
                if isnumeric(SubjectArray(i).subject)
                    z = [ModelInstances.subject] == SubjectArray(i).subject;
                else
                    z = strcmp(SubjectArray(i).subject, {ModelInstances.subject});
                end
                [ModelInstances(z).data] = deal(SubjectArray(i).filename);
            end
        end
        [ModelInstances.data_varname] = deal(p.Results.data_varname);
        [ModelInstances.metadata] = deal(p.Results.metadata);
        [ModelInstances.metadata_varname] = deal(p.Results.metadata_varname);
        bracket_index = 1;
    end

    xyz = cell(numel(SubjectArray),1);
    for i = 1:numel(SubjectArray)
        xyz{i} = SubjectArray(i).getCoords(p.Results.orientation,'xyz','simplify',true);
    end
    switch upper(p.Results.regularization)
        case 'SOSLASSO'
            if isfield(ModelInstances, 'diameter')
                ModelInstances = SOSGroupMake(ModelInstances, xyz);
            elseif isfield(ModelInstances, 'sosgroups')
                ModelInstances = SOSGroupLoad(ModelInstances);
            end
            [ModelInstances.subject] = deal([SubjectArray.subject]);

        case {'LASSO','RIDGE'}
            G = coordGrouping(xyz, 0, 0, 'unitary');
            for ii = 1:numel(ModelInstances)
                ModelInstances(ii).G = G;
            end
    end

    if p.Results.SearchWithHyperband
        n = [p.Results.BRACKETS.n,1];
        r = p.Results.BRACKETS.r;
        while 1
            opts.max_iter = r(bracket_index) * p.Results.IterationsPerHyperband;
            ModelInstances = learn_encoding(ModelInstances, SubjectArray, p.Results.regularization, 'options', opts);
            % Delete low ranked configurations:
            bracket_index = bracket_index + 1;
            ModelInstances = hyperband_pick_top_n(ModelInstances, n(bracket_index));
            if bracket_index < numel(n)
                % save('checkpoint.mat', 'ModelInstances', 'bracket_index');
            else
                if exist('checkpoint.mat', 'file'), delete('checkpoint.mat'); end
                break
            end
        end
    else
        % Grid search
        ModelInstances = learn_encoding(ModelInstances, SubjectArray, p.Results.regularization, 'options', opts);
    end

    cur = 0;
    n = numel(ModelInstances);
    for i = 1:n
        MODEL = ModelInstances(i).Model;
        CONTEXT = rmfield(ModelInstances(i), 'Model');
        if i == 1
            % Preallocate on first pass
            [results,nResultsPerModel] = MODEL.getResults(CONTEXT,SubjectArray,'Initialize',n);
        end
        a = cur + 1;
        b = cur + nResultsPerModel;
        results(a:b) = MODEL.getResults(CONTEXT,SubjectArray);
        cur = b;
    end
    
    %% Save results
    fprintf('Saving stuff...\n');
    rinfo = whos('results');
    switch p.Results.SaveResultsAs
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
end

function [hyperparameters] = verify_setup_MVPA(p)
    % Each regularization requires different lambda configurations. This
    % private function ensures that everything has been properly specified.
    switch upper(p.regularization)
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
            hyperparameters = struct('alpha',alpha,'lambda',lambda,'hyperband',p.SearchWithHyperband);

        case 'LASSO'
            assert(any(isfield(p,{'lambda','lamL1'})), 'Either lambda or lamL1 (synonymns for LASSO regularization) must be defined.');
            if isfield(p,'alpha') && ~isempty(p.alpha)
                warning('Alpha is not relevant for performing Lasso with SOSLasso_logistic. Forcing lamSOS=0, which means that lamL1 will not be scaled up or down. The critical thing for lasso is that all voxels get their own group. This is enforced elsewhere...');
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
                    warning('lamSOS is not relevant for performing Lasso with SOSLasso_logistic. Forcing lamSOS=0 to prevent SOS regularization.');
                    lamSOS = 0;
                else
                    lamSOS = p.lamSOS;
                end
            end
            hyperparameters = struct('lamSOS',lamSOS,'lamL1',lamL1,'hyperband',p.SearchWithHyperband);
            
        case 'RIDGE'
            assert(any(isfield(p,'lamL2')), 'For RIDGE regularization, lamL2 must be defined.');
            if isfield(p,'alpha') && ~isempty(p.alpha)
                warning('Alpha is not relevant for performing Lasso with SOSLasso_logistic. Forcing lamSOS=0, which means that lamL1 will not be scaled up or down. The critical thing for lasso is that all voxels get their own group. This is enforced elsewhere...');

            end
            if isfield(p,'lamSOS') && ~isempty(p.lamSOS)
                if p.lamSOS ~= 0
                    warning('lamSOS is not relevant for performing RIDGE with SOSLasso_logistic. Forcing lamSOS=0 to prevent SOS regularization.');
                end
            end
            lamSOS = 0;
            lamL1 = 0;
            hyperparameters = struct('lamSOS',lamSOS,'lamL1',lamL1,'lamL2',p.lamL2,'hyperband',p.SearchWithHyperband);

        case 'SOSLASSO'
            [b,c] = soslassoParameterExistenceCheck(p);
            assert(b, 'Invalid parameterization for SOS Lasso.');
            if any(c.id == [1,3])
                [lamSOS,lamL1] = ratio2independent(p.alpha,p.lambda);
            else
                lamSOS = p.lamSOS;
                lamL1 = p.lamL1;
            end
            if any(c.id == [3,4])
                if all(size(p.diameter) == [1,3])
                    diameter = {{p.diameter}};
                else
                    diameter = p.diameter;
                end
                if all(size(p.overlap) == [1,3])
                    overlap = {{p.overlap}};
                else
                    overlap = p.overlap;
                end
                hyperparameters = struct('lamSOS',lamSOS,'lamL1',lamL1,'diameter',diameter,'shape',p.shape,'overlap',overlap,'hyperband',p.SearchWithHyperband);
            else
                hyperparameters = struct('lamSOS',lamSOS,'lamL1',lamL1,'sosgroups',p.sosgroups,'hyperband',p.SearchWithHyperband);
            end
            
    end
end

function [hyperparameters] = verify_setup_RSA(p)
    % Each regularization requires different lambda configurations. This private
    % function ensures that everything has been properly specified.
    switch upper(p.regularization)
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
            hyperparameters = struct('alpha',1,'lambda',lam,'hyperband',p.SearchWithHyperband);

        case 'L1L2'
            if isfield(p,'lambda1') && ~isempty(p.lambda1)
                warning('Group Lasso does not use the lambda1 parameter. It is being ignored.');
            end
            if isfield(p,'LambdaSeq') && ~isempty(p.LambdaSeq) && all(~stcmp(p.LambdaSeq,'none'))
                warning('Group Lasso does not use a lambda sequence. Setting to ''none''.');
            end
            assert(~isempty(p.lambda), 'Group Lasso requires lambda.');
            lam    = p.lambda;
            lam1   = NaN;
            lamSeq = 'none';
            hyperparameters = struct('lambda',lam,'lambda1',lam1,'lambdaSeq',lamSeq,'hyperband',p.SearchWithHyperband);

        case {'GROWL','GROWL2'}
            assert(isfield(p,'lambda') && ~isempty(p.lambda) && any(~isnan(p.lambda)), 'grOWL requires lambda.');
            assert(isfield(p,'lambda1') && ~isempty(p.lambda1) && any(~isnan(p.lambda1)), 'grOWL requires lambda1.');
            assert(isfield(p,'LambdaSeq') && ~isempty(p.LambdaSeq), 'A LambdaSeq type (linear or exponential) must be set when using grOWL*.');
            lam    = p.lambda;
            lam1   = p.lambda1;
            lamSeq = p.LambdaSeq;
            hyperparameters = struct('lambda',lam,'lambda1',lam1,'lambdaSeq',lamSeq,'hyperband',p.SearchWithHyperband);

    end
end

function b = islogicallike(x)
    b = islogical(x) || any(x == [1,0]);
end

function b = isintegerlike(x)
    b = mod(x,1) == 0;
end

function b = isMatOrJSON(x)
    b = any(strcmpi(x, {'mat','json'}));
end
