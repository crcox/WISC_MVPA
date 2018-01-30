function WholeBrain_MVPA(varargin)
    p = inputParser;
    p.KeepUnmatched = false;
    % ----------------------Set parameters-----------------------------------------------
    % Model definition
    addParameter(p , 'regularization'   , []        , @ischar        );
    addParameter(p , 'bias'             , false     , @islogicallike );
    addParameter(p , 'lambda'           , []                         );
    addParameter(p , 'alpha'            , []        , @isnumeric     );
    addParameter(p , 'diameter'         , []        , @isnumeric     );
    addParameter(p , 'overlap'          , []        , @isnumeric     );
    addParameter(p , 'shape'            , []        , @ischar        );
    addParameter(p , 'debias'           , false      , @islogicallike );
    addParameter(p , 'AdlasOpts'        , struct()  , @isstruct      );
    % Target definition
    addParameter(p , 'target'           , []        , @ischar        );
    addParameter(p , 'target_type'      , []        , @ischar        );
    % Data definition
    addParameter(p , 'filters'          , []        , @ischarlike    );
    addParameter(p , 'data'             , []        , @ischarlike    );
    addParameter(p , 'data_var'         , 'X'       , @ischar        );
    addParameter(p , 'metadata'         , []        , @ischar        );
    addParameter(p , 'metadata_var'     , 'metadata', @ischar        );
    addParameter(p , 'finalholdout'     , 0         , @isintegerlike );
    addParameter(p , 'cvscheme'         , []        , @isintegerlike );
    addParameter(p , 'cvholdout'        , []        , @isnumeric     );
    addParameter(p , 'orientation'      , []        , @ischar        );
    % Normalization
    addParameter(p , 'normalize'        , false                      );
    % Permutation
    addParameter(p , 'PermutationTest'  , false   ,   @islogicallike );
    addParameter(p , 'PermutationMethod', 'manual'  , @ischar        );
    addParameter(p , 'PermutationIndex' , 'PERMUTATION_STRUCT.mat', @ischar);
    addParameter(p , 'RandomSeed'       , 0                          );
	addParameter(p , 'RestrictPermutationByCV', false, @islogicallike);
    % Hyperband (an alternative to grid search for hyperparameter selection)
    addParameter(p , 'HYPERBAND'        , []                         );
    addParameter(p , 'BRACKETS'         , []                         );
    % Output control
    addParameter(p , 'SmallFootprint'   , false     , @islogicallike );
    addParameter(p , 'SaveResultsAs'    , 'mat'     , @isMatOrJSONOrCSV);
    addParameter(p , 'subject_id_fmt'   , '%d'      , @ischar        );
    % Debugging (none of which currently work...)
    addParameter(p , 'debug'            , false     , @islogicallike );
    addParameter(p , 'SanityCheckData'  , []        , @ischar        );
    % --- searchlight specific --- %
    addParameter(p , 'searchlight'      , false     , @islogicallike );
    addParameter(p , 'slclassifier'     , 'gnb_searchmight',  @ischar);
    addParameter(p , 'slradius'         , []        , @isnumeric     );
    addParameter(p , 'slTestToUse'      , 'accuracyOneSided_analytical', @ischar);
    addParameter(p , 'slpermutations'   , 0         , @isscalar      );
    % Parallel only influences the GLMNET operations, and should only be used
    % when running locally. DO NOT USE ON CONDOR.
    addParameter(p , 'PARALLEL'         , false   ,   @islogicallike );
    % Parameters in this section are unused in the analysis, may exist in
    % the parameter file because other progams use them.
    addParameter(p , 'environment'      , 'condor'  , @ischar        );
    addParameter(p , 'COPY'             , []                         );
    addParameter(p , 'URLS'             , []                         );
    addParameter(p , 'executable'       , []                         );
    addParameter(p , 'wrapper'          , []                         );

    if nargin > 0
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

    % private function.
    required = {'regularization','data','metadata','cvscheme','cvholdout','finalholdout'};
    assertRequiredParameters(p.Results,required);

    DEBUG            = p.Results.debug;
    SmallFootprint   = p.Results.SmallFootprint;
    regularization   = p.Results.regularization;
    debias           = p.Results.debias;
    normalize        = p.Results.normalize;
    BIAS             = p.Results.bias;
    filter_labels    = p.Results.filters;
    target_label     = p.Results.target;
    target_type      = p.Results.target_type;
    datafile         = p.Results.data;
    data_var         = p.Results.data_var;
    cvscheme         = p.Results.cvscheme;
    cvholdout        = p.Results.cvholdout;
    orientation      = p.Results.orientation;
    diameter         = p.Results.diameter;
    overlap          = p.Results.overlap;
    shape            = p.Results.shape;
    slclassifier     = p.Results.slclassifier;
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
    PARALLEL         = p.Results.PARALLEL;
    RandomSeed       = p.Results.RandomSeed;
    PermutationTest  = p.Results.PermutationTest;
    PermutationMethod  = p.Results.PermutationMethod;
    PermutationIndex   = p.Results.PermutationIndex;
    % --- HYPERBAND ---
    HYPERBAND = p.Results.HYPERBAND;
    BRACKETS = p.Results.BRACKETS;

    RestrictPermutationByCV = p.Results.RestrictPermutationByCV;
    SaveResultsAs    = p.Results.SaveResultsAs;
    FMT_subjid       = p.Results.subject_id_fmt;

    if numel(RandomSeed) == 1 && RandomSeed > 0
        rng(RandomSeed);
    end

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
    [X,subjix] = loadData(datafile, data_var, rowfilter, colfilter, metadata, FMT_subjid);
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
        switch PermutationMethod
%                 case 'simple'
%                     if RestrictPermutationByCV
%                         C = permute_target(C, PermutationMethod, cvind);
%                     else
%                         C = permute_target(C, PermutationMethod);
%                     end
            case 'manual'
                load(PermutationIndex, 'PERMUTATION_STRUCT');
                PERMUTATION_INDEX = cell(size(Y));
                for i = 1:numel(Y)
                    %  This is kind of a hack to handle the fact eliminating
                    %  outlying rows and rows belonging to the final holdout
                    %  set will create gaps in the index.
                    [~, ix] = sort(PERMUTATION_STRUCT(i).permutation_index(rowfilter{i}, RandomSeed));
                    [~, permutation_index] = sort(ix);
                    PERMUTATION_INDEX{i} = permutation_index;
                end
            otherwise
                error('crcox:NotImplemented', 'Permutations need to be specified manually.');
        end
    else
        PERMUTATION_INDEX = cell(size(Y));
        for i = 1:numel(Y)
            PERMUTATION_INDEX{i} = (1:size(Y{i}, 1))';
        end
    end
    %% Report whether a bias unit will be used
    fprintf('%-26s', 'Including Bias Unit');
    msg = 'NO';
    if BIAS
        msg = 'YES';
    end
    fprintf(': [%3s]\n', msg);

    % Report whether and how voxels will be standardized
    fprintf('%-26s', 'Normalizing columns of X');
    msg = 'NO';
    if normalize
        msg = normalize;
    end
    fprintf(': [%3s]\n', msg);

    % Final holdout index
    fprintf('%-26s', 'Final holdout index');
    fprintf(': [%3d]\n', finalholdoutInd);

    fprintf('Data loaded and processed.\n');

    %% Plug in the parameters and run
    switch lower(regularization)
        case 'smlr'
            [results,info] = learn_category_encoding(Y, X, regularization, ...
                'lambda'         , lambda         , ...
                'alpha'          , alpha          , ...
                'cvind'          , cvind          , ...
                'cvholdout'      , cvholdout      , ...
                'normalize'      , normalize      , ...
                'bias'           , BIAS           , ...
                'DEBUG'          , DEBUG          , ...
                'SmallFootprint' , SmallFootprint , ...
                'PermutationTest', PermutationTest, ...
                'PermutationMethod', PermutationMethod, ...
                'RestrictPermutationByCV', RestrictPermutationByCV, ...
                'AdlasOpts'      , opts); %#ok<ASGLU>

        case 'lasso_glmnet'
            if isempty(gcp('nocreate')) && PARALLEL
                ppp = parpool('local');
            end

            [results,info] = learn_category_encoding(Y, X, regularization, ...
                'lambda'         , lambda         , ...
                'alpha'          , alpha          , ...
                'cvind'          , cvind          , ...
                'cvholdout'      , cvholdout      , ...
                'normalize'      , normalize      , ...
                'bias'           , BIAS           , ...
                'DEBUG'          , DEBUG          , ...
                'SmallFootprint' , SmallFootprint , ...
                'PARALLEL'       , PARALLEL       , ...
                'permutations'   , PERMUTATION_INDEX, ... % new
                'AdlasOpts'      , opts); %#ok<ASGLU>
            if ~isempty(gcp('nocreate')) && PARALLEL && (exist('ppp', 'var') == 1)
                delete(ppp);
            end
            for ii = 1:numel(results)
%                 [M,ix] = selectbyfield(metadata,'subject',results(ii).subject);
                ix = results(ii).subject;
                M = metadata(ix);
                COORDS = selectbyfield(M.coords, 'orientation', orientation);
                results(ii) = addMaskedCoordinates(results(ii), COORDS, orientation, colfilter{ix});
                results(ii).subject = M.subject;
            end

        case 'iterlasso_glmnet'
            if isempty(gcp('nocreate')) && PARALLEL
                ppp = parpool('local');
            end
            [results,info] = learn_category_encoding(Y, X, regularization, ...
                'lambda'         , lambda         , ...
                'alpha'          , alpha          , ...
                'cvind'          , cvind          , ...
                'cvholdout'      , cvholdout      , ...
                'normalize'      , normalize      , ...
                'bias'           , BIAS           , ...
                'DEBUG'          , DEBUG          , ...
                'SmallFootprint' , SmallFootprint , ...
                'PARALLEL'       , PARALLEL       , ...
                'PermutationTest', PermutationTest, ...
                'PermutationMethod', PermutationMethod, ...
                'RestrictPermutationByCV', RestrictPermutationByCV, ...
                'AdlasOpts'      , opts); %#ok<ASGLU>
            if ~isempty(gcp('nocreate')) && PARALLEL && (exist('ppp', 'var') == 1)
                delete(ppp);
            end
            for ii = 1:numel(results)
                results(ii).subject = metadata(results(ii).subject).subject;
                [M,ix] = selectbyfield(metadata,'subject',results(ii).subject);
                COORDS = selectbyfield(M.coords, 'orientation', orientation);
                results(ii) = addMaskedCoordinates(results(ii), COORDS, orientation, colfilter{ix});
                for iii = 1:numel(results(ii).iterations)
                    results(ii).subject = metadata(results(ii).iterations(iii).subject).subject;
                    results(ii).iterations(iii) = addMaskedCoordinates(results(ii).iterations(iii), COORDS, orientation, colfilter{ix});
                end
            end

        case 'lasso'
            xyz = cell(numel(metadata),1);
            for ii = 1:numel(xyz)
                z = strcmp({metadata(ii).coords.orientation}, orientation);
                xyz{ii} = metadata(ii).coords(z).xyz(colfilter{ii},:);
            end
            % Puts all voxels in one group, which nullifies local and multitask aspects
            % of SOS Lasso, reducing it to lasso. Alpha no longer matters and is set to
            % 0.
            noG = coordGrouping(xyz, inf, 0, 'cube');
            [results,info] = learn_category_encoding(Y, X, regularization, ...
                'groups'         , noG            , ...
                'lambda'         , lambda         , ...
                'alpha'          , 0              , ...
                'cvind'          , cvind          , ...
                'cvholdout'      , cvholdout      , ...
                'normalize'      , normalize      , ...
                'bias'           , BIAS           , ...
                'DEBUG'          , DEBUG          , ...
                'SmallFootprint' , SmallFootprint , ...
                'permutations'   , PERMUTATION_INDEX, ... % new
                'AdlasOpts'      , opts); %#ok<ASGLU>
            
            %% Revise cv indexes
            % Add the final holdout index to all results.
            [results.finalholdout] = deal(finalholdoutInd);
            % Adjust the cvholdout indexes to accomodate the final holdout index.
            %    if isfield(results,'cvholdout') && finalholdoutInd > 0
            %      cvholdout = [results.cvholdout];
            %      z = cvholdout >= finalholdoutInd;
            %      cvholdout(z) = cvholdout(z) + 1;
            %      cvholdout = mat2cell(cvholdout(:),ones(numel(cvholdout),1));
            %      [results.cvholdout] = deal(cvholdout{:});
            %    end
            %% Add extra parameter info
            [results.diameter] = deal(diameter);
            [results.overlap] = deal(overlap);
            [results.shape] = deal(shape);
            for ii = 1:numel(results)
%                 [M,ix] = selectbyfield(metadata,'subject',results(ii).subject);
                ix = results(ii).subject;
                M = metadata(ix);
                COORDS = selectbyfield(M.coords, 'orientation', orientation);
                results(ii) = addMaskedCoordinates(results(ii), COORDS, orientation, colfilter{ix});
                results(ii).subject = M.subject;
            end

        case 'searchlight'
            X = uncell(X);
            Yo = uncell(Y);
            permutations = uncell(PERMUTATION_INDEX);
            nperm = size(permutations,2);
            if min(Yo) == 0
                Yo = Yo + 1;
            end
            results = repmat( ...
                struct('accuracy_map', [], ...
                    'hitrate_map', [], ...
                    'falsealarm', [], ...
                    'pvalue_map', [], ...
                    'subject', [], ...
                    'target', [], ...
                    'RandomSeed', []), ...
                nperm, 1);
            for permix = 1:nperm
                permutation_index = permutations(:,permix);
                Y = Yo(permutation_index);
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
                classifier = slclassifier;
                if strcmp(slTestToUse,'accuracyOneSided_permutation')
                    TestToUseCfg = {'testToUse',slTestToUse,slpermutations};
                else
                    TestToUseCfg = {'testToUse',slTestToUse};
                end
                [am,pm,hm,fm] = computeInformationMap(X,Y,cvind,classifier,'searchlight', ...
                    meta.voxelsToNeighbours,meta.numberOfNeighbours,TestToUseCfg{:});

                results(permix).accuracy_map = am;
                results(permix).hitrate_map = hm;
                results(permix).falsealarm_map = fm;
                results(permix).pvalue_map = pm;
                results(permix).subject = metadata.subject;
                results(permix).target = target_label;
                results(permix).RandomSeed = RandomSeed(permix);
            end

        case 'soslasso'
            xyz = cell(numel(metadata),1);
            for ii = 1:numel(xyz)
                z = strcmp({metadata(ii).coords.orientation}, orientation);
                xyz{ii} = metadata(ii).coords(z).xyz(colfilter{ii},:);
            end
            G = coordGrouping(xyz, diameter, overlap, shape);

            if isempty(HYPERBAND)
                [results,info] = learn_category_encoding(Y, X, regularization, ...
                    'groups'         , G              , ...
                    'lambda'         , lambda         , ...
                    'alpha'          , alpha          , ...
                    'cvind'          , cvind          , ...
                    'cvholdout'      , cvholdout      , ...
                    'normalize'      , normalize      , ...
                    'bias'           , BIAS           , ...
                    'DEBUG'          , DEBUG          , ...
                    'debias'         , debias         , ...
                    'SmallFootprint' , SmallFootprint , ...
                    'permutations'   , PERMUTATION_INDEX, ... % new
                    'AdlasOpts'      , opts); %#ok<ASGLU>

            else
                % NB: Subject and Permutation loop should be handled
                % here for hyperband... this is a ToDo.
%                     [n, r] = hyperband_cfg(HYPERBAND.budget, HYPERBAND.aggressiveness);
%                     s_max = floor((log(HYPERBAND.budget)/log(HYPERBAND.aggressiveness)) + 1);
%                     s = s_max - BRACKETS.s;
                n = BRACKETS.n;
                r = BRACKETS.r;
                SOSLassoInstances = [];
                for i = 1:numel(n)
                    lambda = lambda(1:n(i));
                    opts.max_iter = r(i) * 1000;
                    if ~isempty(SOSLassoInstances)
                        z = ismember([SOSLassoInstances.lambda], lambda);
                        SOSLassoInstances = SOSLassoInstances(z);
                    end
                    [results,SOSLassoInstances] = learn_category_encoding(Y, X, regularization, ...
                        'groups'         , G              , ...
                        'lambda'         , lambda         , ...
                        'alpha'          , alpha          , ...
                        'cvind'          , cvind          , ...
                        'cvholdout'      , cvholdout      , ...
                        'normalize'      , normalize      , ...
                        'bias'           , BIAS           , ...
                        'DEBUG'          , DEBUG          , ...
                        'debias'         , debias         , ...
                        'SmallFootprint' , SmallFootprint , ...
                        'permutations'   , PERMUTATION_INDEX, ... % new
                        'AdlasOpts'      , opts           , ...
                        'SOSLassoInstances' , SOSLassoInstances, ...
                        'hyperband', true);

                    err1 = zeros(n(i), 1);
                    for j = 1:numel(lambda)
                        err1(j) = mean([results([results.lambda] == lambda(j) & [results.alpha] == alpha(j)).err1]);
                    end
                    [~,ix] = sort(err1);
                    lambda = lambda(ix);
                    alpha = alpha(ix);
                end
            end
            %% Revise cv indexes
            % Add the final holdout index to all results.
            [results.finalholdout] = deal(finalholdoutInd);
            % Adjust the cvholdout indexes to accomodate the final holdout index.
            %    if isfield(results,'cvholdout') && finalholdoutInd > 0
            %      cvholdout = [results.cvholdout];
            %      z = cvholdout >= finalholdoutInd;
            %      cvholdout(z) = cvholdout(z) + 1;
            %      cvholdout = mat2cell(cvholdout(:),ones(numel(cvholdout),1));
            %      [results.cvholdout] = deal(cvholdout{:});
            %    end
            %% Add extra parameter info
            [results.diameter] = deal(diameter);
            [results.overlap] = deal(overlap);
            [results.shape] = deal(shape);
            for ii = 1:numel(results)
%                 [M,ix] = selectbyfield(metadata,'subject',results(ii).subject);
                ix = results(ii).subject;
                M = metadata(ix);
                COORDS = selectbyfield(M.coords, 'orientation', orientation);
                results(ii) = addMaskedCoordinates(results(ii), COORDS, orientation, colfilter{ix});
                results(ii).subject = M.subject;
            end
    end
    [results.target] = deal(target_label);
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
            if SmallFootprint
                fields = fieldnames(results);
                fieldIsEmpty = false(1,numel(fields));
                for iField = 1:numel(fields)
                    fn = fields{iField};
                    if all(cellfun(@isempty, {results.(fn)}))
                        fieldIsEmpty(iField) = true;
                    end
                end
                results = rmfield(results,fields(fieldIsEmpty));
            end
            savejson('',results,'FileName','results.json','ForceRootName',false);
        case 'json_testOnly'
            if SmallFootprint
                fields = fieldnames(results);
                fieldIsEmpty = false(1,numel(fields));
                for iField = 1:numel(fields)
                    fn = fields{iField};
                    if all(cellfun(@isempty, {results.(fn)}))
                        fieldIsEmpty(iField) = true;
                    end
                end
                for iResult = 1:numel(results)
                    results(iResult).diff1 = (results(iResult).h1/results(iResult).nt1) - (results(iResult).f1/results(iResult).nd1);
                end
                results = rmfield(results,[fields(fieldIsEmpty)',{'h1','h2','f1','f2','nt1','nt2','nd1','nd2','err2'}]);
            end
            savejson('',results,'FileName','results.json','ForceRootName',false);
        case 'csv'
            % This drops all non-scalar fields, so it is basically small footprint.
            fields = fieldnames(results);
            fieldIsScalar = false(1,numel(fields));
            for iField = 1:numel(fields)
                fn = fields{iField};
                if all(cellfun(@isscalar, {results.(fn)}))
                    fieldIsScalar(iField) = true;
                end
            end
            scalarFields = fields(fieldIsScalar);
            results = rmfield(results,fields(~fieldIsScalar));
            writetable(struct2table(results),'results.csv');
        case 'csv_testOnly'
            % This drops all non-scalar fields, so it is basically small footprint.
            fields = fieldnames(results);
            fieldIsScalar = false(1,numel(fields));
            for iField = 1:numel(fields)
                fn = fields{iField};
                if all(cellfun(@isscalar, {results.(fn)}))
                    fieldIsScalar(iField) = true;
                end
            end
            for iResult = 1:numel(results)
                results(iResult).diff1 = (results(iResult).h1/results(iResult).nt1) - (results(iResult).f1/results(iResult).nd1);
            end
            results = rmfield(results,[fields(~fieldIsScalar)',{'h1','h2','f1','f2','nt1','nt2','nd1','nd2','err2'}]);
            writetable(struct2table(results),'results.csv');
    end
    fprintf('Done!\n');
end

        %% Local functions
function [lambda, alpha] = verifyLambdaSetup(regularization, lambda, alpha)
    % Each regularization requires different lambda configurations. This private
    % function ensures that everything has been properly specified.
    switch lower(regularization)
        case 'smlr'
            if ~isempty(alpha)
                warning('SMLR does not use the alpha parameter. It is being ignored.');
            end
            assert(~isempty(lambda) , 'Lasso requires lambda.');
            alpha  = [];

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
        assert(isfield(params,req), '%s must exist in params structure! Exiting.', req);
        assert(~isempty(params.(req)), '%s must be set. Exiting.', req);
    end
end

function b = isMatOrJSONOrCSV(x)
    b = any(strcmpi(x, {'mat','json','json_testOnly','csv','csv_testOnly'}));
end

function R = addMaskedCoordinates(R, COORDS, ORIENT, COLFILTER)
    C = struct('orientation', ORIENT, 'ind', [], 'ijk', [], 'xyz', []);
    if strcmpi(ORIENT, 'orig')
        if isfield(COORDS, 'ijk')
            ijk = COORDS.ijk(COLFILTER,:);
            C.ijk = ijk(R.Wix,:);
        end
        if isfield(COORDS, 'ind') && ~isempty(COORDS.ind);
            ind = COORDS.ind(COLFILTER);
            C.ind = ind(R.Wix);
        end
    end
    xyz = COORDS.xyz(COLFILTER,:);
    C.xyz = xyz(R.Wix,:);
    R.coords = C;
end
