classdef SOSLasso
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
    addParameter(p , 'verbose' , false                    );
        X
        Y
        alpha
        lambda
        lamSOS
        lam1
        l2 = 0
        Wz % model weights
        maxiter = 100000
        tol = 1e-8
        G
        group_arr
        groups
        num_tasks
        dimension
        t = 1
        iter
        gamma
        gamma_inc = 2
        objective_loss = zeros(100000,1);
    end

    properties ( Access = private, Hidden = true )
        EMPTY = 0;
        Wz_old % model weights
        t_old = 0
    end

    methods
        function obj = SOSLasso(X,Y,alpha,lambda,G,W0,opts)
            if (nargin == 0)
                obj.EMPTY = 1;
                return
            end
            if (nargin <  4), opts = struct(); end
            obj.X = X;
            obj.Y = Y;
            obj.alpha = alpha;
            obj.lambda = lambda;
            obj.G = G;

            % Setup group indexes
            [Gc, ix]  = commonGrouping(G);
            obj.group_arr = group2mat(Gc);
            obj.groups    = group2lab(Gc);
            obj.num_tasks = numel(obj.X);
            obj.dimension = length(obj.groups);

            % initialize (can provide your own initialization)
            if isempty(W0)
                W0 = zeros(obj.dimension, obj.num_tasks);
            else
                for j = 1:length(X);
                    W0{j} = [W0{j}(:); 0];
                    W0{j} = W0{j}(ix(:,j));
                end
                W0 = cell2mat(W0);
            end
            obj.Wz = W0;
            obj.Wz_old = W0;

            % Add dummy unit
            for j = 1:length(X);
                obj.X{j} = [X{j},zeros(size(X{j},1),1)];
                obj.X{j} = X{j}(:,ix(:,j));
                obj.X{j} = X{j}';
            end

            % Set lamL1 and lamSOS
            [obj.lamSOS, obj.lamL1] = ratio2independent(alpha, lambda);
            % Equivalent to:
            %   lamSOS = lambda * alpha;
            %   lamL1  = lambda * (1 - alpha);

            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
            obj.s = RandStream('mt19937ar','Seed',0);
            obj.EMPTY = 0;
        end

        function obj = train(obj, opts)
            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
            obj = SOSLasso_logistic(obj);
        end

        function x = isempty(obj)
            x = all([obj.EMPTY] == 1);
        end

        function W = getW(obj, verbose)
            W = combineOverlappingWeights(obj.Wz,obj.G,'verbose',verbose);
        end
    end
end

function objective_loss = SOSLasso_logistic(obj, verbosity)
    grad_flag = 0;

    while obj.iter < maxiter
        zeta = (obj.t_old - 1) /obj.t;
        Ws = (1 + zeta) * obj.Wz - zeta * obj.Wz_old;

        % compute function value and gradients of the search point
        [grad, Fs ]  = gradVal_eval(Ws);

        % the Armijo Goldstein line search
        while true
            if lamSOS>0 % CRC redefined the "else" block
                Wzp = soslasso_projection(Ws - grad/obj.gamma,lamSOS/obj.gamma,lamL1,group_arr,groups);
            else
%                 Wzp = Ws - grad/gamma;
                Wzp = soslasso_projection(Ws - grad/obj.gamma, 0, lamL1/obj.gamma,group_arr,groups);
            end
%             % CRC moved conditional within the function ( this did not
%             work, because lamL1 needs to be shrunk by gamma I think ... 
%             Wzp = soslasso_projection(Ws - grad/gamma,lamSOS/gamma,lamL1,group_arr,groups);
            Fzp = funVal_eval(Wzp);

            delta_Wzp = Wzp - Ws;
            nrm_delta_Wzp = norm(delta_Wzp, 'fro')^2;
            r_sum = (nrm_delta_Wzp);

            Fzp_gamma = Fs + sum(sum(delta_Wzp.* grad)) + obj.gamma/2 * nrm_delta_Wzp;

            if (r_sum <=1e-20)
                grad_flag=1; % this shows that, the gradient step makes little improvement
                break;
            end

            if (Fzp <= Fzp_gamma)
                break;
            else
                obj.gamma = obj.gamma * obj.gamma_inc;
            end
        end

        obj.Wz_old = obj.Wz;
        obj.Wz = Wzp;

        if (grad_flag)
            break;
        end
%         if lamSOS>0
%             obj(iter) = Fzp + sos_eval(Wz,group_arr,lamSOS,lamL1);
%         else
%             obj(iter) = Fzp;
%         end
        % CRC moved conditional within the function
        obj.objective_loss(obj.iter) = Fzp + sos_eval(obj.Wz,group_arr,lamSOS,lamL1);

        % convergence check.
        if obj.iter>=2
            if (abs( obj(end) - obj(end-1) ) <= tol*obj(end-1))
                break;
            end
        end

        obj.iter = obj.iter + 1;
        obj.t_old = obj.t;
        obj.t = 0.5 * (1 + (1+ 4 * obj.t^2)^0.5);
    end

    % private functions
    function Wshr = soslasso_projection(W,lamSOS,lamL1,group_arr,groups)
        if lamSOS > 0 % CRC added this conditional and what happens if lamSOS == 0 ...
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
        else
            % step 1: perform soft thresholding
            Wshr = sign(W).*max(abs(W) - lamL1,0);
        end
    end

    function [grad_W, funcVal] = gradVal_eval(W)
        grad_W = zeros(dimension, num_tasks);
        lossValVect = zeros (1 , num_tasks);

        for ii = 1:num_tasks
            [ grad_W(:, ii), lossValVect(:, ii)] = unit_grad_eval( W(:, ii), obj.X{ii}, Y{ii});
        end

        grad_W = grad_W + obj.l2 * 2 * W;

        funcVal = sum(lossValVect) + obj.l2*norm(W,'fro')^2;
    end

    function [funcVal] = funVal_eval (W)
        funcVal = 0;

        for ii = 1: num_tasks
            funcVal = funcVal + unit_funcVal_eval( W(:, ii), obj.X{ii}, Y{ii});
        end

        funcVal = funcVal + obj.l2 * norm(W,'fro')^2;
    end


    % SOS regularizer value
    function [regval] = sos_eval(W,group_arr ,lamSOS,lamL1)
        regval = 0;
        [n,~] = size(group_arr);
        Wtemp = [W;zeros(1,num_tasks)];

        if lamSOS > 0 % CRC added this conditional and what happens if lamSOS == 0 ...
            for ii = 1 : n
                w = Wtemp(unique(group_arr(ii,:)), :);
                w = w.^2;
                regval = regval + lamSOS * sqrt(sum(w(:)));
            end
            regval = regval + lamSOS*lamL1*norm(W(:),1);
        else
            regval = regval + lamL1*norm(W(:),1);
        end
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
