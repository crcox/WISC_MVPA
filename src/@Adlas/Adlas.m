classdef Adlas
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        LambdaSequence
        A % Data
        B % Targets
        X % Weights
        trainingFilter
        max_iter = 100000
        fid = 1
        optimIter = 1
        gradIter = 20
        tolInfeas = 1e-6
        tolRelGap = 1e-8
        iter = 0;
        status
        message
        objPrimal
        objDual
        infeas
        trainingError
        testError
    end
    properties ( Access = private, Hidden = true )
        ModelHasBiasUnit
        lambda
        lambda1
        EMPTY = 0;
        AxPrev
        tPrev
        Aprods
        s % RandStream('mt19937ar','Seed',0)
        nitems  % size(A,1)
        nvoxels % size(A,2)
        num_tasks
        r % size(B,2)
        L % Lipschitz constant
        t = 1
        eta = 2
        verbosity = 0
    end

    methods
        function obj = Adlas(A,B,LambdaSequence,trainingFilter,ModelHasBiasUnit,opts)
            if (nargin == 0)
                obj.EMPTY = 1;
                return
            end
            if (nargin <  4), opts = struct(); end
            [obj.nitems,obj.nvoxels] = size(A);
            obj.r = size(B,2);
            obj.s = RandStream('mt19937ar','Seed',0);
            obj.L = 1;
            obj.A = A;
            obj.B = B;
            obj.num_tasks = 1;
            obj.ModelHasBiasUnit = ModelHasBiasUnit;
            % Ensure that lambda is non-increasing
            if ((length(LambdaSequence) > 1) && any(LambdaSequence(2:end) > LambdaSequence(1:end-1)))
                error('Lambda must be non-increasing.');
            end
            if (LambdaSequence(end) < 0)
                error('Lambda must be nonnegative');
            elseif (LambdaSequence(1) == 0)
                error('Lambda must have at least one nonnegative entry.');
            end
            obj.LambdaSequence = LambdaSequence;
            obj.lambda = opts.lambda;
            obj.lambda1 = opts.lambda1;
            if isempty(trainingFilter)
                obj.trainingFilter = true(obj.nitems, 1);
            elseif numel(trainingFilter) ~= size(A,1);
                error('The trainingFilter must have as many elements as there are targets (i.e., examples in the dataset).');
            else
                obj.trainingFilter = trainingFilter;
            end
            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end

            if ~isfield(opts, 'xInit') || (isempty(opts.xInit))
                obj.X = zeros(obj.nvoxels,obj.r);
            end
            obj.EMPTY = 0;
        end

        function obj = train(obj, opts)
            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
            obj = Adlas1(obj);
        end

        function obj = test(obj)
            x = obj.X;                        % Weights
            a = obj.A(~obj.trainingFilter,:); % Data
            b = obj.B(~obj.trainingFilter,:); % Targets
            obj.testError = nrsa_loss(b,a*x);
            a = obj.A(obj.trainingFilter,:);  % Data
            b = obj.B(obj.trainingFilter,:);  % Targets
            obj.trainingError = nrsa_loss(b,a*x);
        end
        function x = isempty(obj)
            x = all([obj.EMPTY] == 1);
        end

        function [W,nzvox,nvox] = getWeights(obj, varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'subjects', 1:numel(obj.A), @isnumeric); % not used
            addParameter(p, 'forceCell', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'dropBias', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'verbose'  , false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});
            W = obj.X;
            if obj.ModelHasBiasUnit && p.Results.dropBias
                W = W(1:end,:);
            end
            nzvox = nnz(any(W,2));
            nvox = size(W,1);
            if numel(p.Results.subjects) == 1 && ~p.Results.forceCell;
                W = W{1};
            end
        end

        function y = getTarget(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'subjects', 1:obj.num_tasks, @isnumeric);
            addParameter(p, 'subset', 'all', @ischar);
            addParameter(p, 'forceCell', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            y = obj.B;
            switch lower(p.Results.subset)
                case 'all'
                    y = y;
                case {'test','testing','testset','testingset'}
                    y = y(~obj.trainingFilter,:);
                case {'train','training','trainset','trainingset'}
                    y = y(obj.trainingFilter,:);
            end
            if numel(p.Results.subjects) == 1 && ~p.Results.forceCell;
                y = y{1};
            end
        end

        function x = getData(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'subjects', 1:obj.num_tasks, @isnumeric);
            addParameter(p, 'subset', 'all', @ischar);
            addParameter(p, 'forceCell', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            x = obj.A(p.Results.subjects);
            switch lower(p.Results.subset)
                case 'all'
%                     x = x;
                case {'test','testing','testset','testingset'}
                    x = x(~obj.trainingFilter,:);
                case {'train','training','trainset','trainingset'}
                    x = x(obj.trainingFilter,:);
            end
            if numel(p.Results.subjects) == 1 && ~p.Results.forceCell;
                x = x{1};
            end
        end

        function [x,y] = getTrainingData(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'subjects', 1:obj.num_tasks, @isnumeric);
            addParameter(p, 'forceCell', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'transposeX', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            x = obj.A;
            y = obj.B;
            if p.Results.transposeX
                x = obj.X(obj.trainingFilter,:)';
            else
                x = obj.X(obj.trainingFilter,:);
            end
            y = obj.Y(obj.trainingFilter,:);
            if numel(p.Results.subjects) == 1 && ~p.Results.forceCell;
                x = x{1};
                y = y{1};
            end
        end

        function [X,Y] = getTestingData(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'subjects', 1:obj.num_tasks, @isnumeric);
            addParameter(p, 'forceCell', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'transposeX', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            if p.Results.transposeX
                X = obj.X(~obj.trainingFilter,:)';
            else
                X = obj.X(~obj.trainingFilter,:);
            end
            Y = obj.Y(~obj.trainingFilter,:);
            if numel(p.Results.subjects) == 1 && ~p.Results.forceCell;
                X = X{1};
                Y = Y{1};
            end
        end

        function Yz = getPrediction(obj, varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'subjects', 1:obj.num_tasks, @isnumeric);
            addParameter(p, 'subset', 'all', @ischar);
            addParameter(p, 'forceCell', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            w = obj.getWeights(p.Results.subjects); % subjects aren't used yet.
            x = obj.getData(p.Results.subjects,'subset',p.Results.subset);
            Yz = cell(size(x));
            Yz = x * w;
            if numel(p.Results.subjects) == 1 && ~p.Results.forceCell;
                Yz = Yz{1};
            end
        end

        function [r,n] = getResults(obj, varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'modelcontext', struct(), @isstruct);
            addOptional(p, 'metadata', struct(), @isstruct);
            addOptional(p, 'subjects', 1:obj.num_tasks, @isnumeric);
            addParameter(p, 'Initialize', 0, @isnumeric);
            addParameter(p, 'SmallFootprint', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            Uz = obj.getWeights(p.Results.subjects,'dropBias',true,'forceCell',true); % subjects aren't used yet
            nz_rows = cellfun(@(x) any(x,2), Uz, 'UniformOutput', 0);
            nzvox = cellfun(@(x) nnz(x), nz_rows, 'UniformOutput', 0);
            nvox = cellfun(@(x) numel(x), nz_rows, 'UniformOutput', 0);
            Yz = obj.getPrediction('subset','all');
            CONTEXT = p.Results.modelcontext;
            META = p.Results.metadata(p.Results.subjects);
            ix = find(any(Uz, 2));
            COORDS = META.coords;
            COORDS_FIELDS = fieldnames(META.coords);
            for j = 1:numel(COORDS_FIELDS)
                cfield = COORDS_FIELDS{j};
                if any(strcmp(cfield, {'ijk','xyz'})) && ~isempty(META.coords.(cfield))
                    COORDS.(cfield) = META.coords.(cfield)(ix,:);
                elseif any(strcmp(cfield, {'ind'})) && ~isempty(META.coords.(cfield))
                    COORDS.(cfield) = META.coords.(cfield)(ix);
                end
            end

            r = struct( ...
                'Uz'               , Uz , ...
                'Cz'               , Cz , ...
                'target_label'     , CONTEXT.target_label , ...
                'target_type'      , CONTEXT.target_type  , ...
                'sim_source'       , CONTEXT.sim_source, ...
                'sim_metric'       , CONTEXT.sim_metric, ...
                'data'             , CONTEXT.data(:) , ...
                'data_varname'     , CONTEXT.data_varname(:) , ...
                'metadata'         , CONTEXT.metadata(:) , ...
                'metadata_varname' , CONTEXT.metadata_varname(:) , ...
                'subject'          , {META.subject}' , ...
                'cvholdout'        , CONTEXT.cvholdout , ...
                'finalholdout'     , CONTEXT.finalholdout , ...
                'regularization'   , CONTEXT.regularization , ...
                'lambda'           , obj.lambda , ...
                'lambda1'          , obj.lambda1 , ...
                'bias'             , obj.ModelHasBiasUnit , ...
                'normalize_wrt'    , CONTEXT.normalize_wrt , ...
                'normalize_data'   , CONTEXT.normalize_data , ...
                'normalize_target' , CONTEXT.normalize_target , ...
                'nz_rows'          , nz_rows , ...
                'nzvox'            , nzvox , ...
                'nvox'             , nvox , ...
                'coords'           , COORDS(:) , ...
                'err1'             , num2cell(obj.testError(:)) , ...
                'err2'             , num2cell(obj.trainingError(:)) , ...
                'RandomSeed'       , CONTEXT.RandomSeed , ...
                'iter'             , obj.iter );

            n = numel(r);
            if p.Results.Initialize > 0
                r = repmat(structfun(@(x) [], r(1), 'UniformOutput', false),p.Results.Initialize * n, 1);
            end
        end

        function disp(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'outputcontrol', 'bycv', @ischar);
            addOptional(p, 'cvindex', 0, @isnumeric);
            addOptional(p, 'subjects', 1:obj.num_tasks, @isnumeric);
            parse(p, obj, varargin{:});

            cvix = p.Results.cvindex;
            subj = p.Results.subjects;
            switch lower(p.Results.outputcontrol)
                case 'header'
                    fprintf('%8s%8s%8s%8s%8s%8s%8s%8s%8s%16s\n', 'subj','cvindex','lambda','lambda1','err1','err2','nzvox','nvox','iter','status');
                case 'bysubject'
                    [~, nzvox, nvox] = obj.getWeights(p.Results.subjects,'forceCell',true,'dropBias',true);
                    fprintf('%8d%8d%8.2f%8.2f%8.2f%8.2f%8d%8d%8d%16s\n', ...
                        subj,cvix,obj.lambda,obj.lambda1,obj.testError,obj.trainingError,nzvox,nvox,obj.iter,obj.message);

                case 'bycv'
                    fprintf('%8d%8d%8.2f%8.2f%8.2f%8.2f%8d%8d%8d%16s\n', ...
                        0,cvix,obj.lambda,obj.lambda1,mean(obj.testError),mean(obj.trainingError),mean(nzvox),mean(nvox),obj.iter,obj.message);
            end

        end

    end
end

function obj = Adlas1(obj, verbosity)
    % Constants for exit status
    STATUS_RUNNING    = 0;
    STATUS_OPTIMAL    = 1;
    STATUS_ITERATIONS = 2;
    STATUS_ALLZERO    = 3;
    STATUS_MSG = {'Optimal','Iteration limit reached','All weights set to zero'};

    if nargin < 2
        verbosity = 0;
    end
    % Initialize parameters
    X       = obj.X; % Weights
    A       = obj.A(obj.trainingFilter,:); % Data
    B       = obj.B(obj.trainingFilter,:); % Targets
    Ax      = A * X;
    Y       = X;
    status  = STATUS_RUNNING;

    tolInfeas = obj.tolInfeas;
    tolRelGap = obj.tolRelGap;

    % Deal with Lasso case
    modeLasso = (numel(obj.LambdaSequence) == 1);
    if (modeLasso)
        proxFunction = @(v1,v2) proxL1L2(v1,v2);
    else
        proxFunction = @(v1,v2) proxSortedL1L2(v1,v2);
    end

    if (verbosity > 0)
        fprintf(fid,'%5s  %9s %9s  %9s  %9s\n','Iter','||r||_F','Gap','Infeas.','Rel. gap');
    end

    % -------------------------------------------------------------
    % Main loop
    % -------------------------------------------------------------
    while (true)

        % Compute the gradient at f(y)
        if (mod(obj.iter,obj.gradIter) == 0) % Includes first iterations
            r = A*Y - B;
            g = A'*(A*Y-B);
            f = trace(r'*r) / 2;
        else
            r = (Ax + ((obj.tPrev - 1) / obj.t) * (Ax - obj.AxPrev)) - B;
            g = A'*(A*Y-B);
            f = trace(r'*r) / 2;
        end

        % Increment iteration count
        obj.iter = obj.iter + 1;

        % Check optimality conditions
        if ((mod(obj.iter,obj.optimIter) == 0))
            % Compute 'dual', check infeasibility and gap
            if (modeLasso)
                gs = sqrt(sum(g.^2,2));
                ys = sqrt(sum(Y.^2,2));

                infeas = max(norm(gs,inf)-obj.LambdaSequence,0);

                objPrimal = f + obj.LambdaSequence*norm(ys,1);
                objDual   = -f - trace(r'*B);
            else
                gs     = sort(sqrt(sum(g.^2,2)),'descend');
                ys     = sort(sqrt(sum(Y.^2,2)),'descend');
                infeas = max(max(cumsum(gs-obj.LambdaSequence)),0);

                % Compute primal and dual objective
                objPrimal =  f + obj.LambdaSequence'*ys;
                objDual  = -f - trace(r'*B);
            end

            % Format string
            if (verbosity > 0)
                str = sprintf(' %9.2e  %9.2e  %9.2e',objPrimal - objDual, infeas/obj.LambdaSequence(1), abs(objPrimal - objDual) / max(1,objPrimal));
            end

            % Check primal-dual gap
            if ((abs(objPrimal - objDual)/max(1,objPrimal) < tolRelGap)  && ...
                    (infeas < tolInfeas * obj.LambdaSequence(1)) )
                status = STATUS_OPTIMAL;
            end

        else
            str = '';
        end

        if (verbosity > 0)
            if ((verbosity == 2) || ...
                    ((verbosity == 1) && (mod(obj.iter,obj.optimIter) == 0)))
                fprintf(fid,'%5d  %9.2e%s\n', obj.iter,f,str);
            end
        end

        % Stopping criteria
        if (status == 0)
            if (obj.iter >= obj.max_iter)
                status = STATUS_ITERATIONS;
            end
        end

        if (status ~= 0)
            if verbosity > 0
                fprintf(fid,'Exiting with status %d -- %s\n', status, STATUS_MSG{status});
            end
            break;
        end

        % Keep copies of previous values
        obj.AxPrev = Ax;
        xPrev  = X;
        fPrev  = f;
        obj.tPrev  = obj.t;

        % Lipschitz search
        while (obj.L < inf)
            % Compute prox mapping
            X = proxFunction(Y - (1/obj.L)*g, obj.LambdaSequence/obj.L);
            d = X - Y;

            Ax = A*X;%A1*vec(X);
            r = Ax-B;
            f = trace(r'*r)/2;
            q = fPrev + trace(d'*g) + (obj.L/2)*trace(d'*d);

            obj.Aprods = obj.Aprods + 1;

            if (q >= f*(1-1e-12))
                break;
            else
                obj.L = obj.L * obj.eta;
            end
        end

        % Update
        obj.t = (1 + sqrt(1 + 4*obj.t^2)) / 2;
        Y = X + ((obj.tPrev - 1) / obj.t) * (X - xPrev);

        % Check if all weights are set to zero
        if all(Y(:)==0) && obj.iter > 100
            status = STATUS_ALLZERO;
            break
        end
    end

    % Set solution
    obj.X = Y;
    obj.objPrimal = objPrimal;
    obj.objDual   = objDual;
    obj.infeas    = infeas;
    if all(Y(:)==0)
        obj.status = STATUS_ALLZERO;
    else
        obj.status = status;
    end
    obj.message = STATUS_MSG{status};
    obj.Aprods  = obj.Aprods + ceil(obj.iter / obj.gradIter);
end

function x = proxL1L2(Y,lambda)
    tmp = Y;
    r = size(Y,2);
    xtmp = tmp./(repmat(sqrt(sum(tmp.^2,2))+realmin,1,r));
    x = xtmp.*repmat(max(sqrt(sum(tmp.^2,2))-lambda,0),1,r);
end
