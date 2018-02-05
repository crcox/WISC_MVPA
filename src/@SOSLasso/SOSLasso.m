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
        ModelHasBiasUnit
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
        function obj = SOSLasso(X,Y,lamSOS,lamL1,G,trainingFilter,ModelHasBiasUnit,opts)
            if (nargin == 0)
                obj.EMPTY = 1;
                return
            end
            if (nargin < 4), opts = struct(); end
            if iscell(X), obj.X = X; else obj.X = {X}; end
            if iscell(Y), obj.Y = Y; else obj.Y = {Y}; end
            obj.lamSOS = lamSOS;
            obj.lamL1 =lamL1;
            obj.G = G;
            obj.ModelHasBiasUnit = ModelHasBiasUnit;

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
            obj.testError = zeros(1, numel(obj.X));
            obj.trainingError = zeros(1, numel(obj.X));
            for i = 1:numel(obj.X)
                w = obj.W(:,i);
                [x,y] = obj.getTrainingData(i);
                obj.testError(i) = classifier_error(y,x*w);
                [x,y] = obj.getTestingData(i);
                obj.trainingError(i) = classifier_error(y,x*w);
            end
        end

        function x = isempty(obj)
            x = all([obj.EMPTY] == 1);
        end

        function W = getWeights(obj, varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'subjects', 1:numel(obj.X), @isnumeric);
            addParameter(p, 'forceCell', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'dropBias', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'verbose'  , false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            Wall = combineOverlappingWeights(obj.W,obj.G,'verbose',p.Results.verbose);
            W = Wall(p.Results.subjects);
            if obj.ModelHasBiasUnit && p.Results.dropBias
                W = cellfun(@(w) w(1:end,:), W, 'UniformOutput', 0);
            end
            if numel(p.Results.subjects) == 1 && ~p.Results.forceCell;
                W = W{1};
            end
        end

        function r = getResults(obj, varargin)
            addOptional(p, 'subjects', 1:numel(obj.X), @isnumeric);
            addOptional(p, 'metadata', struct(), @isstruct);
            addParameter(p, 'Initialize', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'SmallFootprint', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});
            r = struct( ...
                'Wz'               , [] , ...
                'Yz'               , [] , ...
                'target_label'     , [] , ...
                'target_type'      , [] , ...
                'sim_source'       , [] , ...
                'sim_metric'       , [] , ...
                'data'             , [] , ...
                'data_varname'     , [] , ...
                'metadata'         , [] , ...
                'metadata_varname' , [] , ...
                'subject'          , [] , ...
                'cvholdout'        , [] , ...
                'finalholdout'     , [] , ...
                'regularization'   , [] , ...
                'lamSOS'           , [] , ...
                'lamL1'            , [] , ...
                'alpha'            , [] , ...
                'lambda'           , [] , ...
                'bias'             , [] , ...
                'normalizewrt'     , [] , ...
                'normalize_data'   , [] , ...
                'normalize_target' , [] , ...
                'nz_rows'          , [] , ...
                'nzvox'            , [] , ...
                'nvox'             , [] , ...
                'coords'           , [] , ...
                'nt1'              , [] , ...
                'nt2'              , [] , ...
                'nd1'              , [] , ...
                'nd2'              , [] , ...
                'h1'               , [] , ...
                'h2'               , [] , ...
                'f1'               , [] , ...
                'f2'               , [] , ...
                'err1'             , [] , ...
                'err2'             , [] , ...
                'iter'             , [] );

                Wz = obj.getWeights('dropBias',true,'forceCell',true);
                if ~isempty(p.Results.metadata)
                    if numel(p.Results.metadata) > numel(p.Results.subjects)
                        META = p.Results.metadata(p.Results.subjects);
                    end
                end
                    COORDS = cell(size(Wz));
                    COORDS_FIELDS = fieldnames(p.Results.coords);
                    for i = 1:numel(Wz)
                        ix = find(any(Wz, 2));
                        COORDS{i} = p.Results.coords(i)
                        for j = 1:numel(COORDS_FIELDS)
                            cfield = COORDS_FIELDS{j};
                            if any(strcmp(cfield, {'ijk','xyz'})) && ~isempty(p.Results.coords.(cfield))
                                coords = p.Results.coords.(cfield)(ix,:);
                            elseif any(strcmp(cfield, {'ind'})) && ~isempty(COORDS.(cfield))
                                results(iResult).coords.(cfield) = COORDS.(cfield)(ix);
                            end
                        end
                    results(iResult).Uz = Uz;
                    results(iResult).Cz = A.Model.A * A.Model.X;
                end
                results(iResult).subject = A.subject;
                results(iResult).bias = A.bias;
                results(iResult).nz_rows = any(Uz,2);
                results(iResult).nzv = nnz(results(iResult).nz_rows);
                results(iResult).nvox = numel(results(iResult).nz_rows);
                results(iResult).cvholdout = A.cvholdout;
                results(iResult).finalholdout = finalholdoutInd;
                results(iResult).lambda = A.lambda;
                results(iResult).lambda1 = A.lambda1;
                results(iResult).LambdaSeq = A.LambdaSeq;
                results(iResult).regularization = A.regularization;
                results(iResult).tau = tau;
                results(iResult).normalize_data = A.normalize_data;
                results(iResult).normalize_target = normalize_target;
                results(iResult).normalizewrt = A.normalizewrt;
                results(iResult).data = p.Results.data;
                results(iResult).data_var = p.Results.data_varname;
                results(iResult).metadata = p.Results.metadata;
                results(iResult).metadata_var = p.Results.metadata_varname;
                results(iResult).target_label = p.Results.target_label;
                results(iResult).target_type = p.Results.target_type;
                results(iResult).sim_source = p.Results.sim_source;
                results(iResult).sim_metric = p.Results.sim_metric;
                results(iResult).err1 = A.Model.testError;
                results(iResult).err2 = A.Model.trainingError;
                results(iResult).iter = A.Model.iter;
                results(iResult).RandomSeed = A.RandomSeed;
        %         results(iResult).RandomSeed = p.Results.PermutationIndex;
            end
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
        
        function [X,Y] = getTrainingData(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'subjects', 1:numel(obj.X), @isnumeric);
            addParameter(p, 'forceCell', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'transposeX', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            X = cell(obj.X);
            Y = cell(obj.Y);
            for i = p.Results.subjects
                if p.Results.transposeX
                    X{i} = obj.X{i}(obj.trainingFilter{i},:)';
                else
                    X{i} = obj.X{i}(obj.trainingFilter{i},:);
                end
                Y{i} = obj.Y{i}(obj.trainingFilter{i});
            end
            if numel(p.Results.subjects) == 1 && ~p.Results.forceCell;
                X = X{1};
                Y = Y{1};
            end
        end
        
        function [X,Y] = getTestingData(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'subjects', 1:numel(obj.X), @isnumeric);
            addParameter(p, 'forceCell', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'transposeX', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});
            
            X = cell(obj.X);
            Y = cell(obj.Y);
            for i = p.Results.subjects
                if p.Results.transposeX
                    X{i} = obj.X{i}(~obj.trainingFilter{i},:)';
                else
                    X{i} = obj.X{i}(~obj.trainingFilter{i},:);
                end
                Y{i} = obj.Y{i}(~obj.trainingFilter{i});
            end
            if numel(p.Results.subjects) == 1 && ~p.Results.forceCell;
                X = X{1};
                Y = Y{1};
            end
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
                                fprintf('%8d%8d%8.2f%8.2f%8.2f%8.2f%8d%8d%8d%16s\n', ...
                                    subj,cvix,obj.lamSOS,obj.lamL1,obj.testError(i),obj.trainingError(i),nzvox,nvox,obj.iter,obj.message);
                            end
                        end

                    otherwise
                        nzvox = nnz(obj.W);
                        nvox = numel(obj.W);
                        fprintf('%8d%8d%8.2f%8.2f%8.2f%8.2f%8d%8d%8d%16s\n', ...
                            subj,cvix,obj.lamSOS,obj.lamL1,obj.testError,obj.trainingError,nzvox,nvox,obj.iter,obj.message);
                end
            end
        end
    end
end

function obj = SOSLasso_logistic(obj)
    grad_flag = 0;
    [X,Y] = obj.getTrainingData('TransposeX', true);
    
%     dimension = length(obj.groups);
%     num_tasks = numel(obj.X);

    while obj.iter < obj.maxiter
        zeta = (obj.t_old - 1) /obj.t;
        Ws = (1 + zeta) * obj.W - zeta * obj.W_old;

        obj.iter = obj.iter + 1;

        % compute function value and gradients of the search point
        [grad, Fs ]  = gradVal_eval(Ws, X, Y, obj.lamL2);

        % the Armijo Goldstein line search
        while true
            if obj.lamSOS>0 % CRC redefined the "else" block
                Wzp = soslasso_projection(Ws - grad/obj.gamma,obj.lamSOS/obj.gamma,obj.lamL1,obj.group_arr,obj.groups);
            else
%                 Wzp = lasso_projection(Ws - grad/obj.gamma, obj.lamL1/obj.gamma);
                Wzp = Ws - grad/gamma;
            end
            Fzp = funVal_eval (Wzp, X, Y, obj.lamL2);

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
            obj.objective_loss(obj.iter) = Fzp + sos_eval(Wzp,obj.group_arr,obj.lamSOS,obj.lamL1);
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
end

%% --- private functions ---
function Wshr = soslasso_projection(W,lamSOS,lamL1,group_arr,groups)
    num_tasks = size(W,2);
    % step 1: perform soft thresholding
    W_soft = sign(W).*max(abs(W) - lamL1*lamSOS,0);

    %step 2: perform group soft thresholding
    W_soft = [W_soft; zeros(1,size(W_soft,2))]; % for the dummy
    Wtemp = sum(W_soft.^2,2); %xtemp is now a vector
    Wtemp = sum(Wtemp(group_arr),2);
    Wtemp = sqrt(Wtemp);
    Wtemp = max(Wtemp - lamSOS,0); % this is the multiplying factor
    Wtemp = Wtemp./(Wtemp + lamSOS);
    Wtemp = Wtemp(groups);
    Wtemp = repmat(Wtemp,1,num_tasks);
    Wshr = W_soft(1:end-1,:).*Wtemp;
end

function Wshr = lasso_projection(W,lamL1)
    % step 1: perform soft thresholding
    Wshr = sign(W).*max(abs(W) - lamL1,0);
end

function [grad_W, funcVal] = gradVal_eval(W,X,Y,lamL2)
    [dimension, num_tasks] = size(W);
    grad_W = zeros(dimension, num_tasks);
    lossValVect = zeros (1 , num_tasks);

    for ii = 1:num_tasks
        [ grad_W(:, ii), lossValVect(:, ii)] = unit_grad_eval( W(:, ii), X{ii}, Y{ii});
    end

    % If the lamL2 parameter is > 0, then the gradient and funcVal are scaled.
    grad_W = grad_W + (lamL2 * 2 * W);
    funcVal = sum(lossValVect) + lamL2*norm(W,'fro')^2;
end

function [funcVal] = funVal_eval (W, X, Y, lamL2)
    num_tasks = size(W,2);
    funcVal = 0;

    for ii = 1: num_tasks
        funcVal = funcVal + unit_funcVal_eval( W(:, ii), X{ii}, Y{ii});
    end

    funcVal = funcVal + lamL2 * norm(W,'fro')^2;
end

% SOS regularizer value
function [regval] = sos_eval(W, group_arr, lamSOS, lamL1)
    num_tasks = size(W, 2);
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

% function b = validateW0(w0)
%     b = false;
%     if iscell(w0)
%         if numel(obj.X) == numel(w0)
%             b = true;
%         end
%     else
%         if isempty(w0)
%             b = true;
%         end
%     end
% end

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
    parse(p, Wc, G, varargin{:});

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
