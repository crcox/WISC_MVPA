function [W, obj, I] = iterlasso_glmnet(Xtrain, Xtest, Ytrain, Ytest, alpha, lambda, varargin)
    p = inputParser;
    %% Parse function inputs
    %                name         default     validation
    addRequired(p  , 'Xtrain'                            );
    addRequired(p  , 'Xtest'                             );
    addRequired(p  , 'Ytrain'                            );
    addRequired(p  , 'Ytest'                             );
    addRequired(p  , 'alpha'                             );
    addRequired(p  , 'lambda'                            );
    addParameter(p , 'cvind'   , []                      );
    addParameter(p , 'bias'    , 0          , @isscalar  );
    addParameter(p , 'maxiter' , 10         , @isscalar  );
    addParameter(p , 'stopcrit', 1          , @isscalar  );
    addParameter(p , 'tol'     , 1e-8       , @isscalar  );
    addParameter(p , 'W0'      , []         , @validateW0);
    addParameter(p , 'verbose' , false                   );
    addParameter(p , 'PARALLEL', false                   );
    addParameter(p , 'SmallFootprint', false             );
    parse(p, Xtrain, Xtest, Ytrain, Ytest, alpha, lambda, varargin{:});

    Xtrain    = p.Results.Xtrain;
    Xtest     = p.Results.Xtest;
    Ytrain    = cellfun(@checkY, p.Results.Ytrain, 'Unif', 0);
    Ytest     = cellfun(@checkY, p.Results.Ytest, 'Unif', 0);
    alpha     = p.Results.alpha;
    % lambda    = p.Results.lambda;
    bias      = p.Results.bias;
    maxiter   = p.Results.maxiter;
    tol       = p.Results.tol;
    W0        = p.Results.W0;
    CVIND     = p.Results.cvind;
    % verbose   = p.Results.verbose;
    PARALLEL  = p.Results.PARALLEL;
    STOP_CRIT = p.Results.stopcrit;
    SMALL     = p.Results.SmallFootprint;

    if ~iscell(Xtrain)
        Xtrain = {Xtrain};
    end
    if ~iscell(Ytrain)
        Ytrain = {Ytrain};
    end
    if ~iscell(Xtest)
        Xtest = {Xtest};
    end
    if ~iscell(Ytest)
        Ytest = {Ytest};
    end
    if ~iscell(CVIND)
        CVIND = {CVIND};
    end

    W = cell(numel(Xtrain),1);
    I = cell(numel(Xtrain),1);
    for iSubj = 1:numel(Xtrain)
        iterations = struct();
        xtrain = Xtrain{iSubj};
        xtest = Xtest{iSubj};

        total_vox = size(xtest,2);

        ytrain = Ytrain{iSubj};
        ytest = Ytest{iSubj};

        cinds = unique([ytrain(:);ytest(:)]);
        cinds = cinds(:)'; % force to row vec
        m = numel(cinds);
        if m > 2;
            isMultinomial = 1;
            performanceMetric = 'class';
            modelType = 'multinomial';
            ytrain_m = bsxfun(@eq,ytrain(:),cinds);
            ytest_m = bsxfun(@eq,ytest(:),cinds);
        else
            isMultinomial = 0;
            performanceMetric = 'difference';
            modelType = 'binomial';
            ytrain_m = ytrain;
            ytest_m = ytest;
            m = 1;
        end

        used = false(size(xtrain,2),m);
        unused = true(1,size(xtrain,2));

        cvind = CVIND{iSubj};
        cvset = unique(cvind);
        nfold = numel(cvset);
        if cvset(1) > 1
            cvind = cvind - (cvset(1)-1);
            cvset = cvset - (cvset(1)-1);
        end
        adj = find(diff(cvset)>1);
        if ~isempty(adj);
            for a = adj
                cvind(cvind>a) = cvind(cvind>a) - 1;
            end
        end

        % This bit is not used when the objective is to maximize the
        % difference between hits and false alarms, which has a natural
        % chance level of zero.
        if isMultinomial
            tmp = zeros(1,m);
            b = tabulate(cvind);
            for i = 1:m
                a = tabulate(cvind(ytrain==cinds(i)));
                tmp(i) = mean(a(:,2)) ./ mean(b(:,2));
            end
            chance = mean(tmp);
        else
            a = tabulate(cvind(ytrain==max(cinds)));
            b = tabulate(cvind);
            chance = mean(a(:,2)) / mean(b(:,2)); 
        end

        iter = 0;
        nsCounter = 0;
        while 1
            iter = iter + 1;
            opts = glmnetSet(struct('intr',bias,'thresh', tol, 'weights', W0, 'alpha', alpha, 'lambda', []));
            objcv = cvglmnet(xtrain,ytrain,modelType,opts,performanceMetric,nfold,cvind', PARALLEL, 1);
            lambda_min = objcv.lambda_min;

            ci.upper = objcv.cvm - 2*(objcv.cvsd/sqrt(nfold-1));
            h(iter) = ci.upper(objcv.lambda == lambda_min) > 0;
            if isnan(h(iter)) || h(iter) == 0
                nsCounter = nsCounter + 1;
                if nsCounter >= STOP_CRIT;
                    break
                end
            else
                nsCounter = 0;
            end
            if iter > maxiter
                break
            end

            %% Compute single lasso model for the iteration at lambda_min
            opts = glmnetSet(struct('intr',bias,'thresh', tol, 'weights', W0, 'alpha', alpha, 'lambda', lambda_min));
            obj = glmnet(xtrain,ytrain,modelType,opts);
            if isMultinomial
                [~,ix] = max(glmnetPredict(obj,xtest),[],2);
                yztest_m = bsxfun(@eq, ix, 1:m);

                [~,ix] = max(glmnetPredict(obj,xtrain),[],2);
                yztrain_m = bsxfun(@eq, ix, 1:m);
            else
                yztest_m = glmnetPredict(obj,xtest)>0;
                yztrain_m = glmnetPredict(obj,xtrain)>0;
            end

            %% Compute iteration results
            if iscell(obj.beta)
                wz = cell2mat(obj.beta);
            else
                wz = obj.beta;
            end
            nz_rows = any(wz,2);
            tmp = false(total_vox,1);
            tmp(unused) = nz_rows;
            ix = find(tmp);
            nv = size(wz,1);
            wnz = nnz(nz_rows);
            wz = wz(nz_rows,:);

            h1   = sum((ytest_m > 0) & (yztest_m > 0));
            f1   = sum(~(ytest_m > 0) & (yztest_m > 0));
            err1 = sum((ytest_m > 0) ~= (yztest_m > 0));
            nt1  = sum(ytest_m > 0);
            nd1  = sum(ytest_m <= 0);

            h2   = sum((ytrain_m > 0) & (yztrain_m > 0));
            f2   = sum(~(ytrain_m > 0) & (yztrain_m > 0));
            err2 = sum((ytrain_m > 0) ~= (yztrain_m > 0));
            nt2  = sum(ytrain_m > 0);
            nd2  = sum(ytrain_m <= 0);

            if isMultinomial
                [~,yztest] = max(yztest_m,[],2);
                [~,yztrain] = max(yztrain_m,[],2);
                confusion1 = confusionmat(double(ytest),yztest);
                confusion2 = confusionmat(double(ytrain),yztrain);
            else
                [confusion1,confusion2] = deal([]);
            end

            %% Store iteration results
            iterations(iter).subject = uint32(iSubj);
            if ~SMALL
                iterations(iter).Wz = wz;
                iterations(iter).Wix = uint32(ix);
                iterations(iter).Yz = []; % because train and test are passed in separately
            end

            iterations(iter).Wnz = uint32(wnz);
            iterations(iter).nvox = uint32(nv);
            iterations(iter).coords = [];
            iterations(iter).alpha = alpha;
            iterations(iter).lambda = lambda_min;
            iterations(iter).nt1  = uint16(nt1);
            iterations(iter).nt2  = uint16(nt2);
            iterations(iter).nd1  = uint16(nd1);
            iterations(iter).nd2  = uint16(nd2);
            iterations(iter).h1   = uint16(h1);
            iterations(iter).h2   = uint16(h2);
            iterations(iter).f1   = uint16(f1);
            iterations(iter).f2   = uint16(f2);
            iterations(iter).err1 = uint16(err1);
            iterations(iter).err2 = uint16(err2);
            iterations(iter).confusion1 = confusion1;
            iterations(iter).confusion2 = confusion2;

            %% Prepare for next iteration
            % Log used voxels, separate for each classification.
            if iscell(obj.beta)
                used(unused,1:m) = cell2mat(obj.beta) ~= 0;
            else
                used(unused,1:m) = obj.beta ~= 0;
            end

            % Update the unused filter
            unused = ~any(used,2);

            % Reduce X
            xtrain = Xtrain{iSubj}(:,unused);
            xtest = Xtest{iSubj}(:,unused);
        end

        %% Fit final model with ridge regression
        W{iSubj} = zeros(total_vox,m);
        if any(used(:))
            xtrain = Xtrain{iSubj}(:,any(used,2));
            %             xtest = Xtest{iSubj}(:,any(used,2));
            opts = glmnetSet(struct('intr',bias,'thresh', tol, 'weights', W0, 'alpha', 0, 'lambda', []));
            objcv = cvglmnet(xtrain,ytrain,modelType,opts,performanceMetric,nfold,cvind',PARALLEL);
            lambda_min = objcv.lambda_min;

            opts = glmnetSet(struct('intr',bias,'thresh', tol, 'weights', W0, 'alpha', 0, 'lambda', lambda_min));
            obj = glmnet(xtrain,ytrain,modelType,opts);
            if isMultinomial
                W{iSubj}(any(used,2),:) = cell2mat(obj.beta);
            else
                W{iSubj}(any(used,2),:) = obj.beta;
            end
        end
        I{iSubj} = iterations;
    end
end

function b = validateW0(w0)
    b = false;
    if iscell(w0)
        if numel(X) == numel(w0)
            b = true;
        end
    else
        if isempty(w0)
            b = true;
        end
    end
end

function y = checkY(y)
    if islogical(y)
        y = sign(y-0.5);
    end
end

% function dp = evaluateModelFit(yz,y)
%   cinds = unique(y(:)'); % force to row vec
%   m = numel(cinds);
%   if m > 2
%     y1 = bsxfun(@eq,y(:),cinds);
%   else
%     y1 = y > 0;
%   end
%   yz1 = yz > 0;
%
%   hit1 = sum(yz1 & y1);
%   fpos1 = sum(yz1 & ~y1);
%   nt1 = sum(y1);
%   nd1 = sum(~y1);
%   dp = (hit1./nt1) - (fpos1./nd1);
% end
