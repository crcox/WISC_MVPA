function ModelInstances = learn_similarity_encoding(ModelInstances, Y, X, regularization, varargin)
% TO DO:
% Currently, the function can only handle a single subject. It will also
% only support HYPERBAND for a single hyperparameter, meaning that it can
% only do group lasso. This is because the function was originally written
% to facilitate grid search, and with HYPERBAND the hyperparameters should
% be be "crossed" in the same way. With HYPERBAND, configurations can be
% crossed with subject and cvholdout, and probably permutation index.

    p = inputParser();
    addRequired(p  , 'ModelInstances');
    addRequired(p  , 'Y');
    addRequired(p  , 'X');
    addRequired(p  , 'regularization');
    addParameter(p , 'cvind'                   , []       );
    addParameter(p , 'permutations'            , []       );
    addParameter(p , 'AdlasOpts'               , struct() );
    addParameter(p , 'Verbose'                 , true     );
    parse(p, ModelInstances, Y, X, regularization, varargin{:});

    cvind          = p.Results.cvind;
    permutations   = p.Results.permutations;
    options        = p.Results.AdlasOpts;
    VERBOSE        = p.Results.Verbose;

    Xorig = X;
    Yorig = Y;

%     fprintf('Normalize with respect to: %s\n', strrep(normalize_wrt,'_',' '));
%     fprintf('- Normalizing the data: %s ...\n', normalize_data);
%     fprintf('- Normalizing the target dimensions: %s ...\n', normalize_target);

    for i = 1:numel(ModelInstances)
        subix = ModelInstances(i).subject;
        cvix = ModelInstances(i).cvholdout;
        regularization = ModelInstances(i).regularization;
        normalize_data = ModelInstances(i).normalize_data;
        normalize_target = ModelInstances(i).normalize_target;
        normalize_wrt = ModelInstances(i).normalize_wrt;
        BIAS = ModelInstances(i).bias;

        if iscell(X)
            train_set = cell(numel(X),1);
            for j = 1:numel(X)
                P = selectbyfield(permutations,'subject', subix(j), 'RandomSeed', ModelInstances(i).RandomSeed);
                X{j} = Xorig{j};
                Y{j} = Yorig{j}(P.index,:);
                train_set{j}  = cvind{j} ~= cvix; % CHECK THIS
                switch normalize_wrt
                    case 'all_examples'
                        X{j} = normalize_columns(X{j}, normalize_data);
                        Y{j} = normalize_columns(Y{j}, normalize_target);

                    case 'training_set'
                        X{j} = normalize_columns(X{j}, normalize_data, train_set{j});
                        Y{j} = normalize_columns(Y{j}, normalize_target, train_set{j});
                end
                if BIAS
                    X{j} = [X{j}, ones(size(X,1),1)]; %#ok<AGROW> It's not actually growing, see line 86.
                end
            end
        end

        switch upper(regularization)
            case 'L1L2'
                options.lambda = ModelInstances(i).lambda;
                options.lambda1 = ModelInstances(i).lambda1;
                lamseq = options.lambda;
                
            case {'GROWL','GROWL2'}
                % There is no real distinction between GROWL and GROWL2 anymore, but for continuity I'll make GROWL2 map to this anyway.
                [~,d] = size(X{subix});
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

            otherwise
                error('%s is not an implemented regularization. check spelling', regularization);
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
            switch upper(regularization)
                case {'LASSO','SOSLASSO'}
                    ModelInstances(i).Model = SOSLasso(X,Y,lamSOS,lamL1,G,train_set,BIAS,options);
                case {'L1L2','GROWL','GROWL2'}
                    % LambdaSeq must be a column vector
                    ModelInstances(i).Model = Adlas(X{subix}, Y{subix}, lamseq(:), train_set{subix}, BIAS, options);
            end
            ModelInstances(i).Model = ModelInstances(i).Model.train(options);
        elseif ModelInstances(i).Model.status == 2
            ModelInstances(i).Model = ModelInstances(i).Model.train(options);
        else
        % Do nothing
        end
        ModelInstances(i).Model = ModelInstances(i).Model.test();
        if i == 1, disp(ModelInstances(i).Model, 'header'); end
        disp(ModelInstances(i).Model, 'bysubject', ModelInstances(i).cvholdout, ModelInstances(i).subject);
%         err1 = ModelInstances(i).Model.testError;
%         err2 = ModelInstances(i).Model.trainingError;
%         Unz = cellfun(@(x) nnz(any(x,2)), ModelInstances(i).Model.getWeights());
%         
%         nv = size(ModelInstances(i).Model.X, 1);
%         if isempty(lamL1)
%             lamL1 = nan;
%         end
% 
%         if VERBOSE
%             fprintf('%8d%8.2f%8.2f%8.2f%8.2f%8d%8d%8d%16s\n', ...
%                 cvix,lamSOS,lamL1,mean(err1),mean(err2),mean(Unz),mean(nv),ModelInstances(i).Model.iter,'soslasso');
%         end
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

% OLD PERMUTATION ALGORITHM
% -------------------------
%     PERMUTATION_INDEXES = cell(1, max(cvind));
%     for ic = unique(cvind)'
%         disp(sprintf('Permuting CV %d...', ic));
%         c = C(cvind==ic,:);
%         n = size(c,1);
%         if VERBOSE
%             fprintf('Permuting %d rows of C, independently by its %d columns.\n', n, r);
%             fprintf('First 10 rows of C, before shuffling.\n')
%             disp(c)
%         end
%         permix = randperm(n);
%         C(cvind==ic, :) = c(permix, :);
%         if VERBOSE
%             fprintf('First 10 rows of C, after shuffling.\n')
%             disp(c(permix,:))
%         end
%     end

% COMPARISON AGAINST FULL RANK AND LIMITED RANK SQUARE MATRIX
% ===========================================================
% St = C{subix}*C{subix}';
% if strcmpi(target_type,'similarity')
%     lt1  = logical(tril(true(nnz(test_set)),0));
%     s1   = S{subix}(test_set,test_set);
%     sz1  = Sz(test_set,test_set);
%     st1  = St(test_set,test_set);
%     s1   = s1(lt1);
%     sz1  = sz1(lt1);
%     st1  = st1(lt1);
% 
%     lt2 = logical(tril(true(nnz(train_set)),0));
%     s2  = S{subix}(train_set,train_set);
%     sz2 = Sz(train_set,train_set);
%     st2 = St(train_set,train_set);
%     s2  = s2(lt2);
%     sz2 = sz2(lt2);
%     st2 = st2(lt2);
% end

% COLLECT RESULTS
% ===============
% if ~SMALL
%     results(iii).Uz  = uz;
%     results(iii).Uix = uint32(ix(:)');
%     results(iii).Cz  = Cz;
%     results(iii).Sz  = Sz;
% end
% % Metadata
% results(iii).nzv            = uint32(Unz); % number of nonzero rows
% results(iii).nvox           = uint32(nv); % total number of voxels
% results(iii).subject        = subix;
% results(iii).cvholdout      = icv; % cross validation index
% results(iii).finalholdout   = []; % handled in parent function
% if strcmpi(regularization, 'L1L2_GLMNET')
%     results(iii).lambda1    = info.lambda;
% else
%     results(iii).lambda1    = lam1;
% end
% results(iii).lambda         = lam;
% results(iii).LambdaSeq      = LambdaSeq;
% results(iii).regularization = regularization;
% results(iii).bias           = BIAS;
% results(iii).normalize_data = normalize_data;
% 
% if any(test_set)
%     results(iii).err1 = norm(Ch - Chz,'fro')/norm(Ch,'fro');
% end
% results(iii).err2 = norm(Ct - Ctz,'fro')/norm(Ct,'fro');
% 
% if isempty(lambda)
%     lambda_j = nan;
% else
%     lambda_j = lambda(j);
% end
% if isempty(lambda1)
%     lambda1_k = nan;
%     if strcmpi(regularization, 'L1L2_GLMNET')
%         lambda1_k = info.lambda;
%     end
% else
%     lambda1_k = lambda1(k);
% end
% 
% if strcmpi(regularization, 'L1L2_GLMNET')
%     results(iii).iter = info.npasses;
%     info.iter = info.npasses;
% else
%     results(iii).iter = info.iter;
% end

% Old results computation
% =======================
% % KEY
% % ===
% % Our models assume that S = V*W*V'
% % V  : The fMRI data.
% % W  : A matrix of weights.
% % S  : The true similarity matrix.
% % Sz : The predicted similarity matrix.
% % C  : The square root of S, truncated to the first r columns (low rank assumption)
% % Cz : The predicted C.
% % St : The approximated S, reconstructed from actual C.
% % Uz : The estimated voxel weights, with a r weights per voxel.
% 
% nz_rows = any(Uz,2);
% ix = find(nz_rows);
% nv = size(Uz,1);
% Unz = nnz(nz_rows);
% uz = Uz(ix,:);
% 
% % Prepare to evaluate solutions
% Cz = V(:,ix)*uz;
% Ctz = Cz(train_set,:);
% Chz = Cz(test_set,:);
% 
% if any(test_set)
%     err1 = nrsa_loss(Ch, Chz);
%     %cor1 = nrsa_corr(Ch*Ch', Chz*Chz');
% else
%     err1 = [];
%     %cor1 = [];
% end
% err2 = nrsa_loss(Ct, Ctz);
% %cor2 = nrsa_corr(Ct*Ct', Ctz*Ctz');
