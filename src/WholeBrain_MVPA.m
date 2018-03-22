function WholeBrain_MVPA(varargin)
    p = inputParser;
    p.KeepUnmatched = false;
    % ----------------------Set parameters---------------------------------
    % Model definition
	addParameter(p , 'regularization' , []        , @ischar        );
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
    % Target definition
    addParameter(p , 'target_label'     , [], @ischar );
    addParameter(p , 'target_type'      , [], @ischar );
    addParameter(p , 'sim_source'       , [], @ischar );
    addParameter(p , 'sim_metric'       , [], @ischar );
    addParameter(p , 'tau'              , [], @isnumeric     );
    addParameter(p , 'FiltersToApplyBeforeEmbedding' , [], @(x) ischar(x) || iscellstr(x) );
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
        case {'L1L2','GROWL','GROWL2'};
            HYPERPARAMETERS = verify_setup_RSA(p.Results);
        case {'RIDGE','LASSO','SOSLASSO'};
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
    if p.Results.PermutationTest && strcmpi(p.Results.PermutationMethod,'manual');
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
            error('WholeBrain_MVPA:SubjectSelect','Subject ID %s does not uniquely match a metadata entry.', num2cell(SubjectArray.subject));
        end
        % Pull content from metadata object
        CV = M.cvind(:,p.Results.cvscheme);
        RF = selectbyfield(M.filters, 'label', p.Results.filters, 'dimension', 1);
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
        if ~isempty(p.Results.tau) && p.Results.tau ~= 0
            % If tau is set, generate embeddings from target similarity matrices.
            if isempty(p.Results.FiltersToApplyBeforeEmbedding)
                SubjectArray(i) = SubjectArray(i).generateEmbeddings(p.Results.tau, 'ExtendEmbedding', true);
            else
                SubjectArray(i) = SubjectArray(i).generateEmbeddings(p.Results.tau, 'ExtendEmbedding', true, 'PreFilter', p.Results.FiltersToApplyBeforeEmbedding);
            end
        end
        % Set permutation data. If a permutation test is not being run,
        % RandomSeed == 0 and a dummy permutation index will be encoded
        % that is simply 1:Subject.nTotalExamples.
        if p.Results.PermutationTest && strcmpi(p.Results.PermutationMethod, 'manual')
            P = selectbyfield(permutations, 'subject', SubjectArray.subject);
            SubjectArray(i) = SubjectArray(i).setPermutations(p.Results.PermutationMethod, p.Results.RandomSeed, P.permutation_index);
        else
            SubjectArray(i) = SubjectArray(i).setPermutations('none', p.Results.RandomSeed);
        end
    end

    % Report target infomation
    fprintf('\n');
    fprintf('Target Structure Summary\n');
    fprintf('------------------------\n');
    fprintf('%12s: %s\n', 'target_label', p.Results.target_label);
    fprintf('%12s: %s\n', 'type', p.Results.target_type);
    fprintf('%12s: %s\n', 'sim_source', p.Results.sim_source);
    fprintf('%12s: %s\n', 'sim_metric', p.Results.sim_metric);
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
            'subject'          , SubjectsParameter, ...
            'RandomSeed'       , p.Results.RandomSeed       , ...
            'cvholdout'        , p.Results.cvholdout        , ...
            'finalholdout'     , p.Results.finalholdout     , ...
            'bias'             , p.Results.bias             , ...
            'target_label'     , p.Results.target_label     , ...
            'target_type'      , p.Results.target_type      , ...
            'sim_metric'       , p.Results.sim_metric       , ...
            'sim_source'       , p.Results.sim_source       , ...
            'normalize_data'   , p.Results.normalize_data   , ...
            'normalize_target' , p.Results.normalize_target , ...
            'normalize_wrt'    , p.Results.normalize_wrt    , ...
            'regularization'   , p.Results.regularization   , ...
            'HYPERPARAMETERS'  , HYPERPARAMETERS);
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
            [diameter_u,~,ib] = unique(cat(1,ModelInstances.diameter),'rows');
            [overlap_u,~,ic] = unique(cat(1,ModelInstances.overlap),'rows');
            MIT = unique(table( ...
                ib, ...
                ic, ...
                {ModelInstances.shape}', ...
                'VariableNames', {'diameter','overlap','shape'}));
            MIT.diameter = num2cell(MIT.diameter);
            MIT.overlap = num2cell(MIT.overlap);
            G = cell(size(MIT,1),1);
            for ii = 1:size(MIT,1)
                MIT.diameter{ii} = diameter_u(MIT.diameter{ii},:);
                MIT.overlap{ii} = overlap_u(MIT.overlap{ii},:);
                G{ii} = coordGrouping(xyz, ...
                    MIT.diameter{ii}, ...
                    MIT.overlap{ii}, ...
                    MIT.shape{ii});
            end
            for ii = 1:numel(ModelInstances)
                diameter = ModelInstances(ii).diameter;
                overlap = ModelInstances(ii).overlap;
                shape = ModelInstances(ii).shape;
                if size(MIT,1) > 1
                    z = all(bsxfun(@eq, cell2mat(MIT.diameter), diameter), 2) & ...
                        all(bsxfun(@eq, cell2mat(MIT.overlap), overlap), 2) & ...
                        strcmp(MIT.shape, shape);
                else
                    z = 1;
                end 
                ModelInstances(ii).G = G{z};
                ModelInstances(ii).subject = [SubjectArray.subject];
            end

        case {'LASSO','RIDGE'}
            G = coordGrouping(xyz, 0, 0, 'unitary');
            for ii = 1:numel(ModelInstances)
                ModelInstances(ii).G = G;
            end
    end

    if p.Results.SearchWithHyperband
        n = p.Results.BRACKETS.n;
        r = p.Results.BRACKETS.r;
        while 1
            opts.max_iter = r(bracket_index) * p.Results.IterationsPerHyperband;
            ModelInstances = learn_encoding(ModelInstances, SubjectArray, p.Results.regularization, 'AdlasOpts', opts);
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
        ModelInstances = learn_encoding(ModelInstances, SubjectArray, p.Results.regularization, 'AdlasOpts', opts);
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
                lamSOS = 0;
                lamL1 = 0;
            end
            if isfield(p,'lamSOS') && ~isempty(p.lamSOS)
                if p.lamSOS ~= 0
                    warning('lamSOS is not relevant for performing RIDGE with SOSLasso_logistic. Forcing lamSOS=0 to prevent SOS regularization.');
                    lamSOS = 0;
                else
                    lamSOS = p.lamSOS;
                end
            end
            lamSOS = 0;
            lamL1 = 0;
            hyperparameters = struct('lamSOS',lamSOS,'lamL1',lamL1,'lamL2',p.lamL2,'hyperband',p.SearchWithHyperband);

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
            lam1   = [];
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

%         for i = 1:numel(ModelInstances)
%             if iscell(datafiles) && numel(datafiles) > 1;
%                 ModelInstances(i).data = datafiles(ModelInstances(i).subject);
%             else
%                 ModelInstances(i).data = ascell(datafiles);
%             end
%             if iscell(p.Results.data_varname) && numel(p.Results.data_varname) > 1;
%                 ModelInstances(i).data_varname = p.Results.data_varname(ModelInstances(i).subject);
%             else
%                 ModelInstances(i).data_varname = ascell(p.Results.data_varname);
%             end
%             if iscell(p.Results.metadata) && numel(p.Results.metadata) > 1;
%                 ModelInstances(i).metadata = p.Results.metadata(ModelInstances(i).subject);
%             else
%                 ModelInstances(i).metadata = ascell(p.Results.metadata);
%             end
%             if iscell(p.Results.metadata_varname) && numel(p.Results.metadata_varname) > 1;
%                 ModelInstances(i).metadata_varname = p.Results.metadata_varname(ModelInstances(i).subject);
%             else
%                 ModelInstances(i).metadata_varname = ascell(p.Results.metadata_varname);
%             end
%         end
