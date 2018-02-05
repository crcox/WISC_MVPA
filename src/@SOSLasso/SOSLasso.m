classdef SOSLasso
    properties
        X % Data
        Y % Targets
        W % Weights
        G % Groups
        lamSOS
        lamL1
        lamL2 = 0
        trainingFilter
        testError
        trainingError
        maxiter = 100000
        tol = 1e-8
        num_tasks
        dimension
        iter = 0
        objective_loss = zeros(100000,1);
    end

    properties ( Access = public, Hidden = true )
        EMPTY = 0;
        W_old % model weights
        t = 1
        t_old = 0
        s % random seed
        group_arr
        groups
        gamma = 1
        gamma_inc = 2
    end

    methods
        function obj = SOSLasso(X,Y,lamSOS,lamL1,G,trainingFilter,opts)
            if (nargin == 0)
                obj.EMPTY = 1;
                return
            end
            if (nargin < 4), opts = struct(); end
            obj.X = X;
            obj.Y = Y;
            obj.lamSOS = lamSOS;
            obj.lamL1 =lamL1;
            obj.G = G;

            % Setup group indexes
            [Gc, ix]  = commonGrouping(G);
            obj.group_arr = group2mat(Gc);
            obj.groups    = group2lab(Gc);
            obj.num_tasks = numel(obj.X);
            obj.dimension = length(obj.groups);

            if isempty(trainingFilter)
                obj.trainingFilter = true(size(X,1), 1);
            elseif numel(trainingFilter) ~= size(X,1);
                error('The trainingFilter must have as many elements as there are targets (i.e., examples in the dataset).');
            else
                obj.trainingFilter = trainingFilter;
            end
            % initialize (can provide your own initialization)
            % if isempty(W0)
            W0 = zeros(obj.dimension, obj.num_tasks);
            % else
            %     for j = 1:length(X);
            %         W0{j} = [W0{j}(:); 0];
            %         W0{j} = W0{j}(ix(:,j));
            %     end
            %     W0 = cell2mat(W0);
            % end
            obj.W = W0;
            obj.W_old = W0;

            % Add dummy unit
            for j = 1:length(X);
                obj.X{j} = [X{j},zeros(size(X{j},1),1)];
                obj.X{j} = obj.X{j}(:,ix(:,j));
                obj.X{j} = obj.X{j}';
            end

            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
            obj.setRandomSeed(0);
            obj.EMPTY = 0;
        end

        function obj = setRandomSeed(obj, seed)
            obj.s = RandStream('mt19937ar','Seed',seed);
        end
        function obj = train(obj, opts)
            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
            obj = SOSLasso_logistic(obj);
        end

        function obj = test(obj)
            if iscell(obj.W) || size(obj.W,2) > 1
                obj.testError = zeros(1, numel(obj.W));
                obj.trainingError = zeros(1, numel(obj.W));
                for i = 1:numel(obj.W)
                    w = obj.W{i};       % Weights
                    z = obj.trainingFilter{i};
                    x = obj.X{i}(~z,:); % Data
                    y = obj.Y{i}(~z,:); % Targets
                    obj.testError(i) = classifier_error(y,x*w);
                    x = obj.X{i}(z,:);  % Data
                    w = obj.Y{i}(z,:);  % Targets
                    obj.trainingError(i) = classifier_error(y,x*w);
                end
            else
                w = obj.W;       % Weights
                z = obj.trainingFilter;
                x = obj.X(~z,:); % Data
                y = obj.Y(~z,:); % Targets
                obj.testError = classifier_error(y,x*w);
                x = obj.X(z,:);  % Data
                w = obj.Y(z,:);  % Targets
                obj.trainingError = classifier_error(y,x*w);
            end
        end

        function x = isempty(obj)
            x = all([obj.EMPTY] == 1);
        end

        function W = getW(obj, verbose)
            W = combineOverlappingWeights(obj.W,obj.G,'verbose',verbose);
        end

        function group_arr = getGroupArray(obj)
            group_arr = obj.group_arr;
        end

        function group_vec = getGroupVector(obj)
            group_vec = obj.groups;
        end

        function [alpha,lambda] = getAlphaLambda(obj)
            [alpha,lambda] = independent2ratio(obj.lamSOS, obj.lamL1);
        end

        function disp(obj,varargin)
            if nargin < 2
                header = false;
                cvix = 0;
                subject = 0;
            else
                switch lower(varargin{1})
                    case 'header'
                        fprintf('%8s%8s%8s%8s%8s%8s%8s%8s%16s\n', 'subj','cv','lamSOS','lamL1','err1','err2','nzvox','nvox','iter','status');
                    case 'bysubject'
                        cvix = varargin{2};
                        subj = varargin{3};
                        if iscell(obj.W)
                            nvox = zeros(1, numel(obj.W));
                            nzvox = zeros(1, numel(obj.W));
                            for i = 1:numel(obj.W)
                                nzvox = nnz(obj.W{i});
                                nvox = numel(obj.W{i});
                                testError = obj.testError(i);
                                trainingError = obj.trainingError(i);
                                fprintf('%8d%8d%8.2f%8.2f%8.2f%8.2f%8d%8d%8d%16s\n', ...
                                    subj,cvix,obj.lamSOS,obj.lamL1,testError,trainingError,nzvox,nvox,obj.iter,obj.message);
                            end
                        end

                    otherwise
                        nzvox = nnz(obj.W);
                        nvox = numel(obj.W);
                        testError = obj.testError;
                        trainingError = obj.trainingError;
                        fprintf('%8d%8d%8.2f%8.2f%8.2f%8.2f%8d%8d%8d%16s\n', ...
                            subj,cvix,obj.lamSOS,obj.lamL1,testError,trainingError,nzvox,nvox,obj.iter,obj.message);
                end
            end
        end
    end
end

function obj = SOSLasso_logistic(obj)
    grad_flag = 0;
    dimension = length(obj.groups);
    num_tasks = numel(obj.X);
    X = cell(num_tasks,1);
    Y = cell(num_tasks,1);
    for i = 1:num_tasks
        X{i} = obj.X{i}(obj.trainingFilter{i},:);
        Y{i} = obj.Y{i}(obj.trainingFilter{i},:);
    end

    while obj.iter < obj.maxiter
        zeta = (obj.t_old - 1) /obj.t;
        Ws = (1 + zeta) * obj.W - zeta * obj.W_old;

        obj.iter = obj.iter + 1;

        % compute function value and gradients of the search point
        [grad, Fs ]  = gradVal_eval(Ws);

        % the Armijo Goldstein line search
        while true
            if obj.lamSOS>0 % CRC redefined the "else" block
                Wzp = soslasso_projection(Ws - grad/obj.gamma,obj.lamSOS/obj.gamma,obj.lamL1,obj.group_arr,obj.groups);
            else
%                 Wzp = lasso_projection(Ws - grad/obj.gamma, obj.lamL1/obj.gamma);
                Wzp = Ws - grad/gamma;
            end
            Fzp = funVal_eval(Wzp);

            delta_Wzp = Wzp - Ws;
            nrm_delta_Wzp = norm(delta_Wzp, 'fro')^2;
            r_sum = (nrm_delta_Wzp);

            Fzp_gamma = Fs + sum(sum(delta_Wzp.* grad)) + obj.gamma/2 * nrm_delta_Wzp;

            if (obj.iter>1) && (r_sum <=1e-20)
                grad_flag=1; % this shows that, the gradient step makes little improvement
                break;
            end

            if (Fzp <= Fzp_gamma)
                break;
            else
                obj.gamma = obj.gamma * obj.gamma_inc;
            end
        end

        obj.W_old = obj.W;
        obj.W = Wzp;

        if (grad_flag)
            break;
        end

        if obj.lamSOS>0
            obj.objective_loss(obj.iter) = Fzp + sos_eval(obj.W,obj.group_arr,obj.lamSOS,obj.lamL1);
        else
%             obj.objective_loss(obj.iter) = Fzp + L1_eval(obj.W,obj.lamL1);
            obj(iter) = Fzp;
        end

        % convergence check.
        if obj.iter>=2
            if (abs( obj.objective_loss(end) - obj.objective_loss(end-1) ) <= obj.tol*obj.objective_loss(end-1))
                break;
            end
        end

        obj.t_old = obj.t;
        obj.t = 0.5 * (1 + (1+ 4 * obj.t^2)^0.5);
    end

    % private functions
    function Wshr = soslasso_projection(W,lamSOS,lamL1,group_arr,groups)
        % step 1: perform soft thresholding
        X_soft = sign(W).*max(abs(W) - lamL1*lamSOS,0);

        %step 2: perform group soft thresholding
        X_soft = [X_soft; zeros(1,size(X_soft,2))]; % for the dummy
        Xtemp = sum(X_soft.^2,2); %xtemp is now a vector
        Xtemp = sum(Xtemp(group_arr),2);
        Xtemp = sqrt(Xtemp);
        Xtemp = max(Xtemp - lamSOS,0); % this is the multiplying factor
        Xtemp = Xtemp./(Xtemp + lamSOS);
        Xtemp = Xtemp(groups);
        Xtemp = repmat(Xtemp,1,numel(obj.X));
        Wshr = X_soft(1:end-1,:).*Xtemp;
    end

    function Wshr = lasso_projection(W,lamL1)
        % step 1: perform soft thresholding
        Wshr = sign(W).*max(abs(W) - lamL1,0);
    end

    function [grad_W, funcVal] = gradVal_eval(W)
        grad_W = zeros(dimension, num_tasks);
        lossValVect = zeros (1 , num_tasks);

        for ii = 1:num_tasks
            [ grad_W(:, ii), lossValVect(:, ii)] = unit_grad_eval( W(:, ii), obj.X{ii}, obj.Y{ii});
        end

        % If the lamL2 parameter is > 0, then the gradient and funcVal are scaled.
        grad_W = grad_W + (obj.lamL2 * 2 * W);
        funcVal = sum(lossValVect) + obj.lamL2*norm(W,'fro')^2;
    end

    function [funcVal] = funVal_eval (W)
        funcVal = 0;

        for ii = 1: num_tasks
            funcVal = funcVal + unit_funcVal_eval( W(:, ii), obj.X{ii}, obj.Y{ii});
        end

        funcVal = funcVal + obj.lamL2 * norm(W,'fro')^2;
    end


    % SOS regularizer value
    function [regval] = sos_eval(W,group_arr ,lamSOS,lamL1)
        regval = 0;

        [n,~] = size(group_arr);
        Wtemp = [W;zeros(1,num_tasks)];
        for ii = 1 : n
            w = Wtemp(unique(group_arr(ii,:)), :);
            w = w.^2;
            regval = regval + lamSOS * sqrt(sum(w(:)));
        end
        regval = regval + lamSOS*lamL1*norm(W(:),1);
    end

    function [regval] = L1_eval(W,lamL1)
        regval = lamL1*norm(W(:),1);
    end

    function b = validateW0(w0)
        b = false;
        if iscell(w0)
            if numel(obj.X) == numel(w0)
                b = true;
            end
        else
            if isempty(w0)
                b = true;
            end
        end
    end
end

function [ grad_w, funcVal ] = unit_grad_eval( w, x, y)
    %gradient and logistic evaluation for each task
    m = length(y);
    weight = ones(m, 1)/m;
    weighty = weight.* y;
    aa = -y.*(x'*w);
    bb = max( aa, 0);
    funcVal = weight'* ( log( exp(-bb) +  exp(aa-bb) ) + bb );
    pp = 1./ (1+exp(aa));
    b = -weighty.*(1-pp);
    grad_w = x * b;
end

function [ funcVal ] = unit_funcVal_eval( w, x, y)
    %function value evaluation for each task
    m = length(y);
    weight = ones(m, 1)/m;
    aa = -y.*(x'*w);
    bb = max( aa, 0);
    funcVal = weight'* ( log( exp(-bb) +  exp(aa-bb) ) + bb );
end

function y = checkY(y)
    if islogical(y)
        y = sign(y-0.5);
    end
end

function W = combineOverlappingWeights(Wc, G, varargin)
    p = inputParser;
    addRequired(p, 'Wc');
    addRequired(p, 'G');
    addParameter(p , 'verbose' , true);
    parse(p, Wc,G,varargin{:});

    Wc = p.Results.Wc;
    G = p.Results.G;
    verbose = p.Results.verbose;
    % Create indexes into the replicated weight matrix.
    S = subjectfilters(G);
    % Combine (and, in the process, sort into subject order).
    N = cellmax(G);
    W = cell(1,size(G,2));
    if verbose
        fprintf('%8s%8s%8s\n','subject','nnz','unique');
    end
    for j = 1:size(G,2)
        W{j} = zeros(N(j),1);
        for i = 1:size(G,1);
            g = G{i,j};
            s = S{i,j};
            %     s = group_arr(k,:)+1; % to account for bias weight
            if ~isempty(g) && ~isempty(s)
                W{j}(g) = W{j}(g) + Wc(s,j);
            end
        end
        s = cell2mat(S(:,j));
        w = Wc(s,j);
        t = cell2mat(G(:,j));
        if verbose
            fprintf('% 8d% 8d% 8d\n', j, nnz(w), numel(unique(t(w~=0))));
        end
    end
end

function S = subjectfilters(G)
    n = cellfun('length', G);
    m = max(n,[],2);
    mc = [0; cumsum(m)];
    S = cell(size(G));
    for j = 1:size(G,2);
        for i = 1:size(G,1);
            a = mc(i) + 1;
            b = mc(i) + n(i,j);
            S{i,j} = (a:b)';
        end
    end
end

function N = cellmax(C)
    C(cellfun('isempty',C)) = {uint32(0)}; %handle empty cells
    N = uint32(max(cellfun(@max, C),[],1));
end
