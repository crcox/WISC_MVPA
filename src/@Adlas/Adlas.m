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
        n % size(A,2)
        r % size(B,2)
        s % RandStream('mt19937ar','Seed',0)
        L % Lipschitz constant
        t = 1
        eta = 2
        iter = 0;
        status
        message
        verbosity = 0
        objPrimal
        objDual
        infeas
        Aprods
        trainingError
        testError
    end
    properties ( Access = private, Hidden = true )
        EMPTY = 0;
        tPrev
        AxPrev
    end

    methods
        function obj = Adlas(A,B,LambdaSequence,trainingFilter,opts)
            if (nargin == 0)
                obj.EMPTY = 1;
                return
            end
            if (nargin <  4), opts = struct(); end
            obj.n = size(A,2);
            obj.r = size(B,2);
            obj.s = RandStream('mt19937ar','Seed',0);
            obj.L = 1;
            obj.A = A;
            obj.B = B;
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
            if isempty(trainingFilter)
                obj.trainingFilter = true(obj.n, 1);
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
                obj.X = zeros(obj.n,obj.r);
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
