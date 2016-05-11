function [W, obj] = smlr_mvpatoolbox(X, Y, lambda, varargin)
  p = inputParser;
  %% Parse function inputs
  %                name         default     validation
  addRequired(p  , 'X'                                 );
  addRequired(p  , 'Y'                                 );
  addRequired(p  , 'lambda'                            );
  addParameter(p , 'maxiter' , 1000       , @isscalar  );
  addParameter(p , 'tol'     , 1e-8       , @isscalar  );
  addParameter(p , 'W0'      , []         , @validateW0);
  addParameter(p , 'verbose' , false                   );
  parse(p, X, Y, lambda, varargin{:});

  X         = p.Results.X;
  Y         = cellfun(@checkY, p.Results.Y, 'Unif', 0);
  lambda    = p.Results.lambda;
  maxiter   = p.Results.maxiter;
  tol       = p.Results.tol;
  W0        = p.Results.W0;
  verbose   = p.Results.verbose;

  if ~iscell(X)
    X = {X};
  end
  if ~iscell(Y)
    Y = {Y};
  end

  W = cell(numel(X),1);
  for iSubj = 1:numel(X)
    x = X{iSubj};
    y = Y{iSubj};
    if isnan(lambda)
      [W{iSubj},obj(iSubj).args, obj(iSubj).log_posterior, obj(iSubj).wasted, obj(iSubj).saved] = smlr(x,y,'lambda',0.1,'tol',tol,'max_iter',maxiter,'w_init',W0);
    else
      [W{iSubj},obj(iSubj).args, obj(iSubj).log_posterior, obj(iSubj).wasted, obj(iSubj).saved] = smlr(x,y,'lambda',lambda,'tol',tol,'max_iter',maxiter,'w_init',W0);
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
