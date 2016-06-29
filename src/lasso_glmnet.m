function [W, obj] = lasso_glmnet(X, Y,alpha, lambda, varargin)
  p = inputParser;
  %% Parse function inputs
  %                name         default     validation
  addRequired(p  , 'X'                                 );
  addRequired(p  , 'Y'                                 );
  addRequired(p  , 'alpha'                             );
  addRequired(p  , 'lambda'                            );
  addParameter(p , 'cvind'   , []                      );
  addParameter(p , 'maxiter' , 1000       , @isscalar  );
  addParameter(p , 'tol'     , 1e-8       , @isscalar  );
  addParameter(p , 'W0'      , []         , @validateW0);
  addParameter(p , 'verbose' , false                   );
  parse(p, X, Y, alpha, lambda, varargin{:});

  X         = p.Results.X;
  Y         = cellfun(@checkY, p.Results.Y, 'Unif', 0);
  alpha     = p.Results.alpha;
  lambda    = p.Results.lambda;
  maxiter   = p.Results.maxiter;
  tol       = p.Results.tol;
  W0        = p.Results.W0;
  CVIND     = p.Results.cvind;
  verbose   = p.Results.verbose;

  if ~iscell(X)
    X = {X};
  end
  if ~iscell(Y)
    Y = {Y};
  end
  if ~iscell(CVIND)
    CVIND = {CVIND};
  end

  W = cell(numel(X),1);
  for iSubj = 1:numel(X)
    x = X{iSubj};
    y = Y{iSubj};
    
    cvind = CVIND{iSubj};
    cvset = unique(cvind);
    nfold = numel(cvset);
    
    cinds = unique(y);
    cinds = cinds(:)'; % force to row vec
    m = numel(cinds);
    if m > 2;
      isMultinomial = 1;
      performanceMetric = 'class';
      modelType = 'multinomial';
    else
      isMultinomial = 0;
      performanceMetric = 'class';
      modelType = 'binomial';
    end

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
    if isnan(lambda)
      opts = glmnetSet(struct('intr',0,'thresh', tol, 'weights', W0, 'alpha', alpha, 'lambda', []));
      objcv = cvglmnet(x,y,modelType,opts,performanceMetric,nfold,cvind');
      lambda = objcv.lambda_min;

%      yz   = x * objcv.glmnet_fit.beta;
%      nlam = size(yz,2);
%      [h1,h2,f1,f2,nt1,nt2,nd1,nd2,dd1] = deal(zeros(1,nlam));
%      for iLambda = 1:nlam;
%        h1(iLambda)  = nnz( y(test_set{iSubject})>0  & (yz(test_set{iSubject})>0)  );
%        h2(iLambda)  = nnz( y(train_set{iSubject})>0 & (yz(train_set{iSubject})>0) );
%        f1(iLambda)  = nnz(~y(test_set{iSubject})>0  & (yz(test_set{iSubject})>0)  );
%        f2(iLambda)  = nnz(~y(train_set{iSubject})>0 & (yz(train_set{iSubject})>0) );
%        nt1(iLambda) = nnz(y(test_set{iSubject})>0);
%        nt2(iLambda) = nnz(y(train_set{iSubject})>0);
%        nd1(iLambda) = nnz(y(test_set{iSubject})<=0);
%        nd2(iLambda) = nnz(y(train_set{iSubject})<=0);
%        dd1(iLambda) = (h1(iLambda)/nt1(iLambda)) - (f1(iLambda)/nd1(iLambda));
%      end
%      [~,lamind] = max(dd1);
%      lambda = objcv.lambda(lamind);
    end
    opts = glmnetSet(struct('intr',0,'thresh', tol, 'weights', W0, 'alpha', alpha, 'lambda', lambda));
    obj(iSubj) = glmnet(x,y,modelType,opts);
    if iscell(obj(iSubj).beta)
      W{iSubj} = cell2mat(obj(iSubj).beta);
    else
      W{iSubj} = obj(iSubj).beta;
    end
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