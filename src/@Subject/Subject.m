classdef Subject
    %SUBJECT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        filename
        subject
        label
        data
        cvscheme
        rowfilters
        colfilters
        targets
        coords
        permutations
    end

    properties (Hidden=true)
        BoxCar
        WindowStart
        WindowSize
        nTotalExamples
        nTotalFeatures
    end

    methods
        % This function should be more stringent. For example, check that
        % required fields are set and check that they are set to sane
        % values.
        function obj = Subject(varargin)
            if nargin == 0
                return
            elseif nargin == 1
                s = varargin{1};
            else
                s = struct(varargin{:});
            end
            fn = fieldnames(s);
            for i = 1:numel(fn)
                disp(fn{i})
                obj.(fn{i}) = s.(fn{i});
            end
        end

        function obj = setFilename(obj, filename)
            obj.filename = filename;
            if hasWindowInfo(filename)
                s = extractWindowInfo(filename);
                obj.BoxCar = s.BoxCar;
                obj.WindowStart = s.WindowStart;
                obj.WindowSize = s.WindowSize;
            end

            % Functions private to this method
            function s = extractWindowInfo(x)
                a = strfind(x,'BoxCar');
                try
                    y = sscanf(x(a:end), 'BoxCar/%d/WindowStart/%d/WindowSize/%d');
                catch
                    y = sscanf(x(a:end), 'BoxCar\\%d\\WindowStart\\%d\\WindowSize\\%d');
                end
                s = struct('BoxCar',y(1),'WindowStart',y(2),'WindowSize',y(3));
            end

            function b = hasWindowInfo(x)
                y = {'BoxCar','WindowSize','WindowStart'};
                b = all(ismember(y,strsplit(x, '/'))) || all(ismember(y,strsplit(x, '\')));
            end
        end

        function b = hasWindowInfo(obj)
            b = ~isempty(obj.Boxcar);
        end

        function obj = setLabel(obj, label)
            obj.label = label;
        end

        function obj = setDataFromFilenameAndLabel(obj)
            % Currently assumes that data is .mat format
            tmp = load(obj.filename, obj.label);
            obj = obj.setData(tmp.(obj.label));
        end

        function obj = setData(obj, x)
            obj.data = x;
            [obj.nTotalExamples,obj.nTotalFeatures] = size(x);
        end

        function obj = setTargets(obj, targets)
            if (size(targets, 1) == 1)
                targets = targets';
            end
            obj.targets = targets;
        end

        function obj = setRowFilters(obj, filters)
            obj = obj.setFilters(filters, 'rowfilters');
        end

        function obj = setColFilters(obj, filters)
            obj = obj.setFilters(filters, 'colfilters');
        end

        function obj = setFilters(obj,filters,field)
            for i = 1:numel(filters);
                % Force row vector
                filters(i).filter = filters(i).filter(:)';
            end
            obj.(field) = filters;
        end

        function obj = setFinalHoldoutFilter(obj, finalholdoutindex)
            FHO = structWithFields(fieldnames(obj.rowfilters));
            FHO.label = 'finalholdout';
            FHO.dimension = 1;
            FHO.filter = obj.cvscheme(:)' ~= finalholdoutindex;
            if any(strcmpi('finalholdout',{obj.rowfilters.label}))
                obj.rowfilters = replacebyfield(obj.rowfilters,FHO,'label','finalholdout','dimension',1);
            else
                obj.rowfilters(end+1) = FHO;
            end
            function s = structWithFields(x)
                if isstruct(x)
                    fields = fieldnames(x);
                elseif iscellstr(x)
                    fields = x;
                end
                n = numel(fields);
                c = [fields(:)';repmat({[]},1,n)];
                s = struct(c{:});
            end
        end

        function obj = setSubjectIDFromFilename(obj, pattern)
            [~,fname,~] = fileparts(obj.filename);
            id = sscanf(fname, pattern);
            if isempty(id)
                error('Failed to extract subject id from data filename %s. Exiting...', obj.filename);
            end
            obj.subject = id;
        end

        function obj = setSubjectIDFromString(obj, string, pattern)
            [~,fname,~] = fileparts(string);
            id = sscanf(fname, pattern);
            if isempty(id)
                error('Failed to extract subject id from data string %s. Exiting...', string);
            end
            obj.subject = id;
        end

        function obj = setCVScheme(obj, cvscheme)
            obj.cvscheme = cvscheme;
        end

        function obj = setPermutations(obj, method, index, varargin)
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
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addRequired(p, 'method', @ischar);
            addRequired(p, 'index', @isnumeric);
            addOptional(p, 'permpool', [], @isnumeric);
            addParameter(p, 'extend', false, @(x) islogical(x) || any(asnumeric(x) == [1,0]));
            parse(p, obj, method, index, varargin{:});

            switch method
                case 'manual'
                    P = struct( ...
                        'method', p.Results.method, ...
                        'RandomSeed',p.Results.index, ...
                        'index',p.Results.permpool(:,p.Results.index));
                case 'none'
                    P = struct( ...
                        'method', p.Results.method, ...
                        'RandomSeed',p.Results.index, ...
                        'index',(1:obj.nTotalExamples)');
                otherwise
                    error('Subject:setPermutations:InvalidMethod', 'Permutations need to be specified manually.');
            end


            P.index = extend_permutation_index(P.index, obj.nTotalExamples);
            obj.permutations = P;

            function permutation_index = extend_permutation_index(permutation_index, target_length)
                remainder = rem(target_length, size(permutation_index, 1));
                if remainder > 0
                    error('permutation_index has fewer rows than nTotalExamples , and nTotalExamples is not evenly divisible by number in permutation_index.');
                else
                    repeatntimes = target_length / size(permutation_index,1);
                    if repeatntimes > 1
                        warning('permutation_index has fewer rows than nTotalExamples, and nTotalExamples is evenly divisible by number in permutation_index. Repeating permutation_index %d times to match.', repeatntimes);
                        CC = cell(repeatntimes, 1);
                        for k = 1:repeatntimes
                            CC{k} = permutation_index + (size(permutation_index, 1) * (k-1));
                        end
                        permutation_index = cell2mat(CC);
                    end
                end
            end
        end

        function obj = setCoords(obj, coords)
            obj.coords = coords;
        end

        function rf = getRowFilter(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            parse(p, obj, varargin{:});

            if ~isempty(p.Results.include)
                RF = selectbyfield(obj.rowfilters, 'label', p.Results.include, 'dimension', 1);
            elseif ~isempty(p.Results.exclude)
                labels = {obj.rowfilters.label};
                z = ~ismember(labels, p.Results.exclude);
                RF = selectbyfield(obj.rowfilters, 'label', labels(z), 'dimension', 1);
            else
                RF = obj.rowfilters;
            end
            if numel(RF) == 0
                RF(1).filter = true(1, size(obj.data, 1));
            end
            rf = all(cat(1, RF.filter),1);
        end

        function cf = getColFilter(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            parse(p, obj, varargin{:});

            if ~isempty(p.Results.include)
                CF = selectbyfield(obj.colfilters, 'label', p.Results.include, 'dimension', 2);
            elseif ~isempty(p.Results.exclude)
                labels = {obj.colfilters.label};
                z = ~ismember(labels, p.Results.exclude);
                CF = selectbyfield(obj.colfilters, 'label', labels(z), 'dimension', 2);
            else
                CF = obj.colfilters;
            end
            if numel(CF) == 0
                CF(1).filter = true(1, size(obj.data, 2));
            end
            cf = all(cat(1, CF.filter),1);
        end

        function [rf,cf] = getFilters(obj, varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            parse(p, obj, varargin{:});

            rf = obj.getRowFilter(varargin{:});
            cf = obj.getColFilter(varargin{:});
        end

        function c = getCoords(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addOptional(p, 'orientation', '', @ischar);
            addOptional(p, 'coordsystem', '', @ischar);
            addParameter(p, 'unfiltered', false, @(x) islogical(x) || any(asnumeric(x) == [1,0]));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'simplify', false, @(x) islogical(x) || any(asnumeric(x) == [1,0]));
            parse(p, obj, varargin{:});

            if p.Results.unfiltered
                cf = true(1, size(obj.data,2));
            else
                if ~isempty(p.Results.include);
                    cf = obj.getColFilter('include',p.Results.include);
                elseif ~isempty(p.Results.exclude);
                    cf = obj.getColFilter('exclude',p.Results.exclude);
                else
                    cf = obj.getColFilter();
                end
            end
            if ~isempty(p.Results.orientation)
                x = selectbyfield(obj.coords,'orientation', p.Results.orientation);
            else
                x = obj.coords;
            end
            if ~isempty(p.Results.coordsystem)
                x.(p.Results.coordsystem) = x.(p.Results.coordsystem)(cf, :);
            else
                for i = 1:numel(x)
                    if isfield(x(i), 'ind') && ~isempty(x(i).ind), x(i).ind = x(i).ind(cf,:); end
                    if isfield(x(i), 'ijk') && ~isempty(x(i).ijk), x(i).ijk = x(i).ijk(cf,:); end
                    if isfield(x(i), 'xyz') && ~isempty(x(i).xyz), x(i).xyz = x(i).xyz(cf,:); end
                end
            end
            if p.Results.simplify && ~isempty(p.Results.coordsystem) && (numel(x) == 1)
                c = x(1).(p.Results.coordsystem);
            else
                c = x;
            end
        end

        function x = getCVScheme(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addParameter(p, 'unfiltered', false, @(x) islogical(x) || any(asnumeric(x) == [1,0]));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            parse(p, obj, varargin{:});

            if p.Results.unfiltered
                x = obj.cvscheme;
            else
                if ~isempty(p.Results.include);
                    rf = obj.getRowFilter('include',p.Results.include);
                elseif ~isempty(p.Results.exclude);
                    rf = obj.getRowFilter('exclude',p.Results.exclude);
                else
                    rf = obj.getRowFilter();
                end
                x = obj.cvscheme(rf);
            end
        end

        function t = getTestSet(obj, cvind, varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addRequired(p, 'cvind', @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'unfiltered', false, @(x) islogical(x) || any(asnumeric(x) == [1,0]));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            parse(p, obj, cvind, varargin{:});

            x = obj.getCVScheme( ...
                'unfiltered', p.Results.unfiltered, ...
                'include', p.Results.include, ...
                'exclude', p.Results.exclude);

            t = x == p.Results.cvind;
        end

        function t = getTrainingSet(obj,cvind,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addRequired(p, 'cvind', @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'unfiltered', false, @(x) islogical(x) || any(asnumeric(x) == [1,0]));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            parse(p, obj, cvind, varargin{:});

            t = ~obj.getTestSet( ...
                p.Results.cvind, ...
                'unfiltered', p.Results.unfiltered, ...
                'include', p.Results.include, ...
                'exclude', p.Results.exclude);
        end

        function x = getData(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addParameter(p, 'unfiltered', false, @(x) islogical(x) || any(asnumeric(x) == [1,0]));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            parse(p, obj, varargin{:});

            if p.Results.unfiltered
                x = obj.data;
            else
                if ~isempty(p.Results.include);
                    [rf, cf] = obj.getFilters('include',p.Results.include);
                elseif ~isempty(p.Results.exclude);
                    [rf, cf] = obj.getFilters('exclude',p.Results.exclude);
                else
                    [rf, cf] = obj.getFilters();
                end
                x = obj.data(rf, cf);
            end
        end

        function c = getTargets(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addParameter(p, 'unfiltered', false, @(x) islogical(x) || any(isnumeric(x) == [1,0]));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'simplify', false, @(x) islogical(x) || any(isnumeric(x) == [2,1,0]));
            parse(p, obj, varargin{:});

            if ~p.Results.unfiltered
                if ~isempty(p.Results.include);
                    rf = obj.getRowFilter('include',p.Results.include);
                elseif ~isempty(p.Results.exclude);
                    rf = obj.getRowFilter('exclude',p.Results.exclude);
                else
                    rf = obj.getRowFilter();
                end
            end

            x = obj.targets;
            for i = 1:numel(x)
                if iscell(x(i).type)
                    for j = 1:numel(x(i).type)
                        if p.Results.unfiltered
                            x(i).target{j} = x(i).target{j};
                        elseif strcmpi(x(i).type{j},'similarity')
                            x(i).target{j} = x(i).target{j}(rf,rf);
                        else
                            x(i).target{j} = x(i).target{j}(rf,:);
                        end
                    end
                else
                    if p.Results.unfiltered
                        x(i).target = x(i).target;
                    elseif strcmpi(x(i).type,'similarity')
                        x(i).target = x(i).target(rf,rf);
                    else
                        x(i).target = x(i).target(rf,:);
                    end
                end
            end
            if p.Results.simplify && (numel(x) == 1)
                if iscell(x(1).target)
                    % This is poorly implemented, but simplify can take on
                    % either a logical value or a numeric value of 0, 1,
                    % *or 2*. If 2, then we select from "the history". By
                    % this, I am refering to the fact that similarity
                    % matrices tagged with type 'similarity' might be
                    % converted to embeddings. In these cases, multiple
                    % target values will be stored in a cell array, with
                    % most recent version (the embedding) being the first
                    % element and the prior version (the similarity matrix)
                    % being the second.
                    % The overloading of the simplicity parameter is
                    % further abused by exploiting that true behaves the
                    % same way as 1 in this context. This means if simplify
                    % == true or simplify == 1, the more recent version of
                    % the target (the embedding) will be selected, and if
                    % simplify == 2 the similarity matrix will be selected.
                    c = x(1).target{p.Results.simplify};
                else
                    c = x(1).target;
                end
            else
                c = x;
            end
        end

        function pix = getPermutationIndex(obj, RandomSeed)
            z = obj.permutations.RandomSeed == RandomSeed;
            pix = obj.permutations.index(:,z);
        end

        function tp = getPermutedTargets(obj,RandomSeed,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addRequired(p, 'RandomSeed', @(x) isnumeric(x) || isscalar(x));
            addParameter(p, 'unfiltered', false, @(x) islogical(x) || any(isnumeric(x) == [1,0]));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'simplify', false, @(x) islogical(x) || any(isnumeric(x) == [2,1,0]));
            parse(p, obj, RandomSeed, varargin{:});

            t = obj.getTargets(...
                'unfiltered', true, ...
                'include', p.Results.include, ...
                'exclude', p.Results.exclude, ...
                'simplify', p.Results.simplify);

            pix = obj.getPermutationIndex(p.Results.RandomSeed);
            if p.Results.unfiltered
                tp = t(pix,:);
            else
                rf = obj.getRowFilter();
                rfp = rf(pix);
                tp = t(pix,:);
                tp = tp(rfp,:);
            end
        end

        function obj = generateEmbeddings(obj, tau, varargin)
        % For all targets with type 'similarity', generate a low-rank
        % embedding and update the 'target' field of the targets structure
        % (replacing the item-by-item symetric similarity matrix.
        % PreFilterWith indicates which row filters should be applied to
        % the similarity matrix before computing the embedding.
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addRequired(p, 'tau', @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'ExtendEmbedding', false, @(x) islogical(x) || any(asnumeric(x) == [1,0]));
            addParameter(p, 'PreFilterWith', [], @(x) ischar(x) || iscellstr(x));
            parse(p, obj, tau, varargin{:});

            if strcmpi(obj.targets.type, 'similarity');
                T = obj.getTargets('unfiltered',true,'simplify',false);
                if isempty(p.Results.PreFilterWith)
                    S = obj.getTargets('unfiltered',true,'simplify',true);
                else
                    S = obj.getTargets('include',p.Results.PreFilterWith,'simplify',true);
                end
                [C, r] = sqrt_truncate_r(S, tau);
                fprintf('S decomposed into %d dimensions (tau=%.2f)\n', r, tau);

                if ~isempty(p.Results.PreFilterWith)
                    tmp = C;
                    rf = obj.getRowFilter('include',p.Results.PreFilterWith);
                    C = nan(numel(rf),r);
                    C(rf,:) = tmp;
                end

                if isfield(T,'embedding_subset') && ~isempty(T.embedding_subset)
                    if ~isempty(p.Results.PreFilterWith)
                        warning('Embedding Subsets and Similarity Matrix Prefilters are probably incompatible features! This may be an error in the future.');
                    end
                    C = C(T.embedding_subset,:);
                end

                C = expand(C, obj.nTotalExamples);
                obj.targets.target = {C,S};
                obj.targets.type = {'embedding','similarity'};
            end
            function C = expand(C,n)
                if size(C,1) < n
                    remainder = rem(n, size(C,1));
                    if remainder > 0
                        error('C has fewer rows than X, and number in X is not evenly divisible by number in C.');
                    else
                        repeatntimes = n ./ size(C,1);
                        warning('C has fewer rows than X, and number in X is evenly divisible by number in C. Repeating C %d times to match.', repeatntimes);
                        C = repmat(C, repeatntimes, 1);
                    end
                end
            end
        end
    end
end
