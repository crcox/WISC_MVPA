function ModelInstances = learn_encoding(ModelInstances, SubjectArray, regularization, varargin)
    p = inputParser();
    addRequired(p  , 'ModelInstances');
    addRequired(p  , 'SubjectArray');
    addRequired(p  , 'regularization');
    addParameter(p , 'options' , struct() );
    parse(p, ModelInstances, SubjectArray, regularization, varargin{:});

    for i = 1:numel(ModelInstances)
        S = selectbyfield(SubjectArray,'subject',ModelInstances(i).subject);
        X = cell(numel(S), 1);
        Y = cell(numel(S), 1);
        train_set = cell(numel(X),1);
        for j = 1:numel(S)
            X{j} = S(j).getData();
            Y{j} = S(j).getPermutedTargets(ModelInstances(i).RandomSeed,'simplify',true);
            train_set{j}  = S(j).getTrainingSet(ModelInstances(i).cvholdout);
            switch ModelInstances(i).normalize_wrt
                case 'all_examples'
                    X{j} = normalize_columns(X{j}, ModelInstances(i).normalize_data);
                    Y{j} = normalize_columns(Y{j}, ModelInstances(i).normalize_target);

                case 'training_set'
                    X{j} = normalize_columns(X{j}, ModelInstances(i).normalize_data, train_set{j});
                    Y{j} = normalize_columns(Y{j}, ModelInstances(i).normalize_target, train_set{j});
            end
            if ModelInstances(i).bias
                X{j} = [X{j}, ones(size(X,1),1)];
            end
        end

        bias = ModelInstances(i).bias;
        options = p.Results.options;
        switch upper(ModelInstances(i).regularization)
            case 'L1L2'
                X = X{1};
                Y = Y{1};
                train_set = train_set{1};
                options.lambda = ModelInstances(i).lambda;
                options.lambda1 = NaN;
                lamseq = options.lambda;

            case {'GROWL','GROWL2'}
            % There is no real distinction between GROWL and GROWL2
            % anymore, but for continuity I'll make GROWL2 map to this
            % anyway.
                X = X{1};
                Y = Y{1};
                train_set = train_set{1};
                d = size(X, 2);
                options.lambda = ModelInstances(i).lambda;
                options.lambda1 = ModelInstances(i).lambda1;
                LambdaSeq = ModelInstances(i).lambdaSeq;
                switch lower(LambdaSeq)
                    case 'linear'
                        lamseq = options.lambda1*(d:-1:1)/d + options.lambda;
                    case 'exponential'
                        lamseq = options.lambda*sqrt(2*log((d*ones(1,d))./(1:d)));
                    case 'inf'
                        lamseq = [options.lambda+options.lambda1, repmat(options.lambda,1,d-1)];
                end
                
            case {'LASSO','SOSLASSO'}
                lamSOS = ModelInstances(i).lamSOS;
                lamL1 = ModelInstances(i).lamL1;
                G = ModelInstances(i).G;
                
            case 'RIDGE'
                lamSOS = ModelInstances(i).lamSOS;
                lamL1 = ModelInstances(i).lamL1;
                options.lamL2 = ModelInstances(i).lamL2;
                G = ModelInstances(i).G;

            otherwise
                error('%s is not an implemented regularization. check spelling', ModelInstances(i).regularization);
        end

        % TRAINING CONDITIONS
        % ===================
        % In case of new model:
        % --------------------
        %   1. Initialize model
        %   2. Train model
        %
        % In case of existing model:
        % -------------------------
        %   1. Check that status == 2, which means that the previous round
        %   of training stopped because it hit the iteration limit.
        %   2. If so, train for more iterations.
        %
        % In case of existing model and status == 1 or status == 3:
        % --------------------------------------------------------
        %   Continue without doing anything. Status 1 means optimal
        %   convergence with at least one nonzero weight assigned. Status 3
        %   means the solution is zero-sparse.
        %
        %   If a zero-sparse solution is obtained in a single iteration,
        %   this can prevent some information from ever being logged in the
        %   Adlas structure, which results in an error at run-time when
        %   trying to pick up where things left off.
        %
        %   By not trying to train these models any more (which is fine,
        %   because they have converged on a solution already, anyway),
        %   this error should be avoided.
        if isempty(ModelInstances(i).Model)
            switch upper(ModelInstances(i).regularization)
                case {'RIDGE','LASSO','SOSLASSO'}
                    ModelInstances(i).Model = SOSLasso(lamSOS,lamL1,G,train_set, bias,options);
                case {'L1L2','GROWL','GROWL2'}
                    ModelInstances(i).Model = Adlas(size(X),size(Y),lamseq(:), train_set, bias, options);
            end
            ModelInstances(i).Model = ModelInstances(i).Model.train(X,Y,options);
            ModelInstances(i).Model = ModelInstances(i).Model.test(X,Y);
            if i == 1, printresults(ModelInstances(i).Model, 'header'); end
            printresults(ModelInstances(i).Model, 'bysubject', ModelInstances(i).cvholdout, ModelInstances(i).subject);
        elseif ModelInstances(i).Model.status == 2
            ModelInstances(i).Model = ModelInstances(i).Model.train(X,Y,options);
            ModelInstances(i).Model = ModelInstances(i).Model.test(X,Y);
            if i == 1, printresults(ModelInstances(i).Model, 'header'); end
            printresults(ModelInstances(i).Model, 'bysubject', ModelInstances(i).cvholdout, ModelInstances(i).subject);
        else
            % Do nothing
        end
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
        warning('There are %d constant-valued voxels. These voxels are not normalized.', sum(~z));
    end
end
