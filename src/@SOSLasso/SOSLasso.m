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
        message = '';
    end

    methods
        function obj = SOSLasso(X,Y,lamSOS,lamL1,G,trainingFilter,ModelHasBiasUnit,opts)
            if (nargin == 0)
                obj.EMPTY = 1;
                return
            end
            if (nargin < 4), opts = struct(); end
            if iscell(X), obj.X = X(:); else obj.X = {X}; end
            if iscell(Y), obj.Y = Y(:); else obj.Y = {Y}; end
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
                obj.trainingFilter = true(size(obj.X{1},1), 1);
            elseif numel(trainingFilter) ~= size(obj.X,1);
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
            obj.testError = zeros(1, obj.num_tasks);
            obj.trainingError = zeros(1, obj.num_tasks);
            for i = 1:obj.num_tasks
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

        function group_arr = getGroupArray(obj)
            group_arr = obj.group_arr;
        end

        function group_vec = getGroupVector(obj)
            group_vec = obj.groups;
        end

        function [alpha,lambda] = getAlphaLambda(obj)
            [alpha,lambda] = independent2ratio(obj.lamSOS, obj.lamL1);
        end

        function [W,nzvox,nvox] = getWeights(obj, varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'subjects', 1:obj.num_tasks, @isnumeric);
            addParameter(p, 'forceCell', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'dropBias', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'combine', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'verbose'  , false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            if p.Results.combine
                Wall = combineOverlappingWeights(obj.W,obj.G,'verbose',p.Results.verbose);
            else
                Wall = mat2cell(obj.W,size(obj.W,1),ones(1,size(obj.W,2)));
            end
            W = Wall(p.Results.subjects);
            if obj.ModelHasBiasUnit && p.Results.dropBias
                W = cellfun(@(w) w(1:end,:), W, 'UniformOutput', 0);
            end
            nzvox = cellfun(@(x) nnz(any(x,2)), W, 'UniformOutput', 1);
            nvox = cellfun(@(x) size(x,1), W, 'UniformOutput', 1);
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

            y = obj.Y(p.Results.subjects);
            switch lower(p.Results.subset)
                case 'all'
                    for i = 1:numel(y), y{i} = y{i}; end
                case {'test','testing','testset','testingset'}
                    for i = 1:numel(y), y{i} = y{i}(~obj.trainingFilter{i}); end
                case {'train','training','trainset','trainingset'}
                    for i = 1:numel(y), y{i} = y{i}(obj.trainingFilter{i}); end
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

            x = obj.X(p.Results.subjects);
            switch lower(p.Results.subset)
                case 'all'
                    for i = 1:numel(x), x{i} = x{i}; end
                case {'test','testing','testset','testingset'}
                    for i = 1:numel(x), x{i} = x{i}(~obj.trainingFilter{i},:); end
                case {'train','training','trainset','trainingset'}
                    for i = 1:numel(x), x{i} = x{i}(obj.trainingFilter{i},:); end
            end
            if numel(p.Results.subjects) == 1 && ~p.Results.forceCell;
                x = x{1};
            end
        end

        function [X,Y] = getTrainingData(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'subjects', 1:obj.num_tasks, @isnumeric);
            addParameter(p, 'forceCell', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'transposeX', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            X = cell(size(obj.X));
            Y = cell(size(obj.Y));
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
            addOptional(p, 'subjects', 1:obj.num_tasks, @isnumeric);
            addParameter(p, 'forceCell', false, @(x) islogical(x) || x==1 || x==0);
            addParameter(p, 'transposeX', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            X = cell(size(obj.X));
            Y = cell(size(obj.Y));
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

        function Yz = getPrediction(obj, varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'subjects', 1:obj.num_tasks, @isnumeric);
            addParameter(p, 'subset', 'all', @ischar);
            addParameter(p, 'forceCell', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            w = obj.getWeights(p.Results.subjects,'forceCell',true,'dropBias',false,'combine',false);
            x = obj.getData(p.Results.subjects,'subset',p.Results.subset);
            Yz = cell(size(x));
            for i = 1:numel(x)
                Yz{i} = x{i} * w{i};
            end
            if numel(p.Results.subjects) == 1 && ~p.Results.forceCell;
                Yz = Yz{1};
            end
        end

        function [truepos,falsepos,poscount,negcount] = getTFPos(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'subjects', 1:obj.num_tasks, @isnumeric);
            addParameter(p, 'subset', 'test', @ischar);
            parse(p, obj, varargin{:});

            x = obj.getTarget('subset',p.Results.subset);
            P = obj.getPrediction('subset',p.Results.subset);

            poscount = zeros(size(x));
            negcount = zeros(size(x));
            truepos = zeros(size(x));
            falsepos = zeros(size(x));
            for i = 1:numel(x)
                pp = P{i} > 0;
%                 pn = P{i} <= 0;
                xp = x{i} > 0;
                xn = x{i} <= 0;

                poscount(i) = nnz(xp);
                negcount(i) = nnz(xn);

                truepos(i) = nnz(xp & pp);
                falsepos(i) = nnz(xn & pp);
            end
        end

        function [r,n] = getResults(obj, varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'modelcontext', struct(), @isstruct);
            addOptional(p, 'metadata', struct(), @(x) isa(x,'Subject'));
            addOptional(p, 'subjects', 1:obj.num_tasks, @isnumeric);
            addParameter(p, 'Initialize', 0, @isnumeric);
            addParameter(p, 'SmallFootprint', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            Wz = obj.getWeights(p.Results.subjects,'dropBias',true,'forceCell',true,'combine',true);
            nz_rows = cellfun(@(x) any(x,2), Wz, 'UniformOutput', 0);
            nzvox = cellfun(@(x) nnz(x), nz_rows, 'UniformOutput', 0);
            nvox = cellfun(@(x) numel(x), nz_rows, 'UniformOutput', 0);
            Yz = obj.getPrediction('subset','all');
            [truepos.testset,falsepos.testset,poscount.testset,negcount.testset] = ...
                obj.getTFPos(p.Results.subjects,'subset','TestSet');
            [truepos.trainingset,falsepos.trainingset,poscount.trainingset,negcount.trainingset] = ...
                obj.getTFPos(p.Results.subjects,'subset','TrainingSet');
            [alpha,lambda] = obj.getAlphaLambda();
            CONTEXT = p.Results.modelcontext;
            META = p.Results.metadata(p.Results.subjects);
            COORDS = cell(size(Wz));
            for i = 1:numel(Wz)
                ix = find(any(Wz{i}, 2));
                COORDS{i} = META(i).coords;
                COORDS_FIELDS = fieldnames(META(i).coords);
                for j = 1:numel(COORDS_FIELDS)
                    cfield = COORDS_FIELDS{j};
                    if any(strcmp(cfield, {'ijk','xyz'})) && ~isempty(META(i).coords.(cfield))
                        COORDS{i}.(cfield) = META(i).coords.(cfield)(ix,:);
                    elseif any(strcmp(cfield, {'ind'})) && ~isempty(META(i).coords.(cfield))
                        COORDS{i}.(cfield) = META(i).coords.(cfield)(ix);
                    end
                end
            end

            r = struct( ...
                'Wz'               , Wz , ...
                'Yz'               , Yz , ...
                'target_label'     , CONTEXT.target_label , ...
                'target_type'      , CONTEXT.target_type  , ...
                'data'             , CONTEXT.data(:) , ...
                'data_varname'     , CONTEXT.data_varname(:) , ...
                'metadata'         , CONTEXT.metadata(:) , ...
                'metadata_varname' , CONTEXT.metadata_varname(:) , ...
                'subject'          , {META.subject}' , ...
                'cvholdout'        , CONTEXT.cvholdout , ...
                'finalholdout'     , CONTEXT.finalholdout , ...
                'regularization'   , CONTEXT.regularization , ...
                'lamSOS'           , obj.lamSOS , ...
                'lamL1'            , obj.lamL1 , ...
                'alpha'            , alpha , ...
                'lambda'           , lambda , ...
                'bias'             , obj.ModelHasBiasUnit , ...
                'normalize_wrt'    , CONTEXT.normalize_wrt , ...
                'normalize_data'   , CONTEXT.normalize_data , ...
                'normalize_target' , CONTEXT.normalize_target , ...
                'nz_rows'          , nz_rows , ...
                'nzvox'            , nzvox , ...
                'nvox'             , nvox , ...
                'coords'           , COORDS(:) , ...
                'nt1'              , num2cell(poscount.testset) , ...
                'nt2'              , num2cell(poscount.trainingset) , ...
                'nd1'              , num2cell(negcount.testset) , ...
                'nd2'              , num2cell(negcount.trainingset) , ...
                'h1'               , num2cell(truepos.testset) , ...
                'h2'               , num2cell(truepos.trainingset) , ...
                'f1'               , num2cell(falsepos.testset) , ...
                'f2'               , num2cell(falsepos.trainingset) , ...
                'err1'             , num2cell(obj.testError(:)) , ...
                'err2'             , num2cell(obj.trainingError(:)) , ...
                'RandomSeed'       , CONTEXT.RandomSeed , ...
                'iter'             , obj.iter );

            n = numel(r);
            if p.Results.Initialize > 0
                r = repmat(structfun(@(x) [], r(1), 'UniformOutput', false),p.Results.Initialize * n, 1);
            end
        end

        function printresults(obj,varargin)
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
                    fprintf('%8s%8s%8s%8s%8s%8s%8s%8s%8s%16s\n', 'subj','cvindex','lamSOS','lamL1','err1','err2','nzvox','nvox','iter','status');
                case 'bysubject'
                    [w, nzvox, nvox] = obj.getWeights(p.Results.subjects,'forceCell',true,'dropBias',true,'combine',true);
                    for i = 1:numel(w)
                        fprintf('%8d%8d%8.2f%8.2f%8.2f%8.2f%8d%8d%8d%16s\n', ...
                            subj(i),cvix,obj.lamSOS,obj.lamL1,obj.testError(i),obj.trainingError(i),nzvox(i),nvox(i),obj.iter,obj.message);
                    end

                case 'bycv'
                    fprintf('%8d%8d%8.2f%8.2f%8.2f%8.2f%8d%8d%8d%16s\n', ...
                        0,cvix,obj.lamSOS,obj.lamL1,mean(obj.testError),mean(obj.trainingError),mean(nzvox),mean(nvox),obj.iter,obj.message);
            end

        end
    end
end

function obj = SOSLasso_logistic(obj)
    grad_flag = 0;
    [X,Y] = obj.getTrainingData('TransposeX', true,'forceCell',true);

%     dimension = length(obj.groups);
%     num_tasks = obj.num_tasks;

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
%         if obj.num_tasks == numel(w0)
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
    W = cell(size(G,2),1);
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
