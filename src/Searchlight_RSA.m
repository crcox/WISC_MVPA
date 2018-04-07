function results = Searchlight_RSA( varargin )
%SEARCHLIGHT_MVPA A wrapper to SearchMight (by Francisco Periera)
    p = Searchlight_MVPA_Parameters();
    p = parse_input_parameters(p, varargin);
    set_global_random_stream_seed(p.Results.RandomSeed)
    
    metadata = load_variable(p.Results.metadata, p.Results.metadata_varname);
    % Load permutation object
    if p.Results.PermutationTest && strcmpi(p.Results.PermutationMethod,'manual');
        permutations = load_variable(p.Results.PermutationIndex, p.Results.perm_varname);
    end
    
    SubjectArray = load_data_as_subjects( ...
        p.Results.data, ...
        p.Results.data_varname, ...
        p.Results.subject_id_fmt);

    SubjectArray = add_metadata_to_subjects( ...
        SubjectArray, ...
        metadata, ...
        p.Results.data, ...
        p.Results.cvscheme, ...
        p.Results.finalholdout, ...
        p.Results.target_label, ...
        p.Results.target_type, ...
        p.Results.sim_source, ...
        p.Results.sim_metric, ...
        p.Results.filters);
    
    for i = 1:numel(SubjectArray)
        if ~isempty(p.Results.tau) && p.Results.tau ~= 0
            % If tau is set, generate embeddings from target similarity matrices.
            if isempty(p.Results.FiltersToApplyBeforeEmbedding)
                SubjectArray(i) = SubjectArray(i).generateEmbeddings(p.Results.tau, 'ExtendEmbedding', true);
            else
                SubjectArray(i) = SubjectArray(i).generateEmbeddings(p.Results.tau, 'ExtendEmbedding', true, 'PreFilter', p.Results.FiltersToApplyBeforeEmbedding);
            end
        end
        if p.Results.PermutationTest && strcmpi(p.Results.PermutationMethod, 'manual')
            P = selectbyfield(permutations, 'subject', SubjectArray(i).subject);
            SubjectArray(i) = SubjectArray(i).setPermutations(p.Results.PermutationMethod, p.Results.RandomSeed, P.permutation_index);
        else
            SubjectArray(i) = SubjectArray(i).setPermutations('none', p.Results.RandomSeed);
        end
    end
    
    report_target_information( ...
        p.Results.target_label, ...
        p.Results.target_type, ...
        p.Results.sim_source, ...
        p.Results.sim_metric)
    
    report_data_information(SubjectArray)
    fprintf('Data loaded and processed.\n');
    
    results = initialize_results_struct(p, SubjectArray);
    for i = 1:numel(SubjectArray)
        SL = run_searchlight_models( ...
            SubjectArray(i), ...
            p.Results.cvholdout, ...
            p.Results.normalize_data, ...
            p.Results.normalize_target, ...
            p.Results.radius, ...
            p.Results.orientation, ...
            p.Results.permutations, ...
            p.Results.SortByCoordsIndex);
        
        if i == 1, printresults(ModelInstances(i).Model, 'header'); end
        printresults(ModelInstances(i).Model, 'bysubject', p.Results.cvholdout, SubjectArray(i).subject);
        results(i) = update_results(results(i),SL);
    end
    save('results.mat', 'results');
end

function p = parse_input_parameters(p, args)
    if isempty(args)
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
        
    elseif isstruct(args{1})
        s = args{1};
        x = [fieldnames(s), struct2cell(s)]';
        parse(p, x{:});
        
    else
        % From command line
        parse(p, args{:});
    end
end

function p = Searchlight_MVPA_Parameters()
    p = inputParser;
    p.KeepUnmatched = false;
    % ----------------------Set parameters---------------------------------
    % Model definition
    addParameter(p , 'classifier' , ''      , @ischar        );
    addParameter(p , 'bias'       , false   , @islogicallike );
    % Searchlight definition
    addParameter(p , 'radius'           , []      , @isnumeric );
    addParameter(p , 'shape'            , []      , @ischar    );
    addParameter(p,  'SortByCoordsIndex', false, @islogicallike);
    % Target definition
    addParameter(p , 'target_label'     , [], @ischar );
    addParameter(p , 'target_type'      , 'category', @ischar );
    addParameter(p , 'sim_source'       , [], @ischar );
    addParameter(p , 'sim_metric'       , [], @ischar );
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
    addParameter(p , 'normalize_wrt'   , 'all_examples', @ischar ); % NEW VARIABLE!
    % Permutation
    addParameter(p , 'RandomSeed'   , 0, @(x) isnumeric(x) && isscalar(x));
    addParameter(p , 'permutations' , 0, @(x) isnumeric(x) && isscalar(x));
    % Output control
    addParameter(p , 'SmallFootprint', false  , @islogicallike );
    addParameter(p , 'SaveResultsAs'  , 'mat' , @isMatOrJSON   );
    addParameter(p , 'subject_id_fmt' , '%d'  , @ischar        );
    % Parameters in this section are unused in the analysis, may exist in
    % the parameter file because other progams use them.
    addParameter(p , 'HYPERBAND'  , [] );
    addParameter(p , 'SearchWithHyperband', false );
    addParameter(p , 'BRACKETS' , [] );
    addParameter(p , 'IterationsPerHyperband', 0, @isnumeric );
    addParameter(p , 'COPY'       , [] );
    addParameter(p , 'URLS'       , [] );
    addParameter(p , 'executable' , [] );
    addParameter(p , 'wrapper'    , [] );
    
    function b = isMatOrJSON(x)
        b = any(strcmpi(x, {'mat','json'}));
    end
    function b = islogicallike(x)
        b = islogical(x) || any(x == [0,1]);
    end
    function b = isintegerlike(x)
        b = mod(x,1) == 0;
    end
end

function set_global_random_stream_seed(seed)
    if ~isempty(seed) && isscalar(seed)
        if seed == 0;
            rng('default')
        else
            rng(seed);
        end
    end
end

function x = load_variable(filename,variable)
    StagingContainer = load(filename, variable);
    x = StagingContainer.(variable);
end

function SubjectArray = load_data_as_subjects(filenames, variable, subject_id_fmt)
    datafiles  = ascell(filenames);
    N = length(datafiles);
    SubjectArray = repmat(Subject, 1, N);
    for i = 1:N
        x = datafiles{i};
        fprintf('Loading %s from  %s...\n', variable, x);
        % Log and process data filename and data label
        SubjectArray(i) = SubjectArray(i).setFilename(datafiles{i});
        SubjectArray(i) = SubjectArray(i).setLabel(variable);
        SubjectArray(i) = SubjectArray(i).setDataFromFilenameAndLabel();
        SubjectArray(i) = SubjectArray(i).setSubjectIDFromFilename(subject_id_fmt);
    end
end

function SubjectArray = add_metadata_to_subjects(SubjectArray, data_filenames, metadata, cvscheme, finalholdout, target_label, target_type, sim_source, sim_metric, filters)
    for i = 1:numel(SubjectArray)
        % Select metadata matching current subject ID
        M = selectbyfield(metadata, 'subject', SubjectArray(i).subject);
        if numel(M) > 1
            error('WholeBrain_MVPA:SubjectSelect','Subject ID %s does not uniquely match a metadata entry.', num2cell(SubjectArray.subject));
        end
        % Pull content from metadata object
        CV = M.cvind(:,cvscheme);
        
        RF = selectbyfield(M.filters, 'label', filters, 'dimension', 1);
        CF = selectbyfield(M.filters, 'label', filters, 'dimension', 2);
        T = selectbyfield(M.targets, ...
            'label', target_label, ...
            'type', target_type, ...
            'sim_source', sim_source, ...
            'sim_metric', sim_metric);

        SubjectArray(i) = SubjectArray(i).setCVScheme(CV);
        SubjectArray(i) = SubjectArray(i).setRowFilters(RF);
        SubjectArray(i) = SubjectArray(i).setColFilters(CF);
        SubjectArray(i) = SubjectArray(i).setFinalHoldoutFilter(finalholdout);
        SubjectArray(i) = SubjectArray(i).setCoords(M.coords);
        SubjectArray(i) = SubjectArray(i).setTargets(T);
    end
end

function report_target_information(target_label, target_type, sim_source, sim_metric)
    % Report target infomation
    fprintf('\n');
    fprintf('Target Structure Summary\n');
    fprintf('------------------------\n');
    fprintf('%12s: %s\n', 'target_label', target_label);
    fprintf('%12s: %s\n', 'type', target_type);
    fprintf('%12s: %s\n', 'sim_source', sim_source);
    fprintf('%12s: %s\n', 'sim_metric', sim_metric);
    fprintf('\n');
end

function report_data_information(SubjectArray)
    fprintf('Data Dimensions\n');
    fprintf('---------------\n');
    fprintf('%16s%16s%16s\n','subject','initial','filtered');
    fprintf('%s\n',repmat('-',1,16*3));
    for i = 1:numel(SubjectArray)
        fprintf('%16s (%6d,%6d) (%6d,%6d)\n', ...
            num2str(SubjectArray(i).subject), ...
            size(SubjectArray(i).getData('unfiltered',true)), ...
            size(SubjectArray(i).getData('unfiltered',false)));
    end
    fprintf('\n');
end

function results = initialize_results_struct(p,SubjectArray)
    results = struct( ...
        'subject', {SubjectArray}, ...
        'target_label', p.Results.target_label, ...
        'condition', p.Results.orientation, ...
        'shape', p.Results.shape, ...
        'radius', p.Results.radius, ...
        'RandomSeed', p.Results.RandomSeed, ...
        'error_map1', [], ...
        'error_map2', []);
end

function r = update_results(r,SL)
    r.error_map1 = SL.error_map1;
    r.error_map2 = SL.error_map2;
end

function SL = run_searchlight_models(S, cvholdout, normalize_data, normalize_target, radius, orientation, RandomSeed, SortByCoordsIndex)
    X = S.getData();
    Y = S.getPermutedTargets(RandomSeed,'simplify',true);
    normalize_wrt = 'all_examples';
    switch normalize_wrt
        case 'all_examples'
            X = normalize_columns(X, normalize_data);
            Y = normalize_columns(Y, normalize_target);

        case 'training_set'
            % Cross validation is done within GLMNET, so we
            % cannot normalize this way.
    end

    coords = selectbyfield(S.coords,'orientation',orientation);
    % SortByCoordsIndex is used to implement manipulations in simulations
    if SortByCoordsIndex
        X = X(:,coords.ind);
        [~,xi] = sort(coords.ind);
    end
    train_set  = S.getTrainingSet(cvholdout);
    xyz = S.getCoords(ModelInstances.orientation, 'xyz', 'simplify', true);
    SL = Searchlight(radius,xyz,S.cvscheme,train_set{1},bias,struct());
    SL.computeInformationMap(X{1},Y{1});
    
    if SortByCoordsIndex
        SL.error_map1 = SL.error_map1(xi);
        SL.error_map2 = SL.error_map2(xi);
    end    
end

function y = normalize_columns(x, method, wrt)
    if nargin < 3
        wrt = true(size(x,1),1);
    end
    % By default, subtract zero from each column.
    mm = zeros(1, size(x,2));
    % By default, divide each column by one.
    ss = ones(1, size(x,2));
    switch lower(method)
        case 'none'
            % Do nothing
        case {'zscore','zscored'}
            mm = mean(x(wrt,:),1);
            ss = std(x(wrt,:),0,1);
        case {'center','centered','centre','centred'}
            mm = mean(x(wrt,:),1);
        case 'stdev'
            ss = std(x(wrt,:),0,1);
        case '2norm'
            mm = mean(x(wrt,:),1);
            ss = norm(x(wrt,:));
        otherwise
            error('Unrecognized normalizaion method "%s"! Exiting...', method);
    end
    % Avoid dividing by zero, when columns have constant value.
    z = ss > 0;
    y(:,z) = bsxfun(@minus,x(:,z), mm(z));
    y(:,z) = bsxfun(@rdivide,y(:,z), ss(z));
    
    if any(~z)
        warning('There are %d constant-valued voxels. These voxels are not normalized.', sum(z));
        if VERBOSE
            fprintf('Constant-valued voxel indexes:\n');
            disp(find(~z));
        end
    end
end