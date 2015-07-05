function [W, obj] = SOG_logistic(X, Y,alpha,lambda,groups,group_arr,varargin)

% function solves
% min \sum_t log(1 + exp(-Yt.*XtW(:,t))) + alpha Omega(W)
%
%
% INPUT
% X         : {n x d} length t cell array
% Y         : {n x 1} length t output
% alpha     : regularizer for the SOSlasso penalty
% lambda    : regularizer for the L1 part of soslasso penalty
% groups    : group assignment (see related function)
% group_arr : group matrix (see related functions)
% optional  :
% l2        : an additional l2 regularizer for smoothing the solution
% maxiter   : max iterations
% tol       : tolerance
%
% OUTPUT
% W: model: d x t

%%%%%%%%%
% parts of this code have been adapted from the MALSAR package
% http://www.public.asu.edu/~jye02/Software/MALSAR/
%%%%%%%%%

  p = inputParser;
  o = inputParser;
  %% Parse function inputs
  %                name         default     validation
  addRequired(p  , 'X'                                );
  addRequired(p  , 'Y'                                );
  addRequired(p  , 'alpha'                , @isscalar );
  addRequired(p  , 'lambda'               , @isscalar );
  addRequired(p  , 'groups'               , @isnumeric);
  addRequired(p  , 'group_arr'            , @isnumeric);
  addParameter(p , 'opts'      , struct() , @isstruct );
  parse(p, X, Y, alpha, lambda, groups, group_arr, varargin{:});

  X         = cellfun(@transpose, p.Results.X, 'Unif', 0);
  Y         = p.Results.Y;
  lambda    = p.Results.lambda;
  alpha     = p.Results.alpha;
  groups    = p.Results.groups;
  group_arr = p.Results.group_arr;
  opts      = p.Results.opts;
  optfields = fieldnames(opts);
  optcell   = [optfields'; struct2cell(opts)'];

  %% Parse optional options structure
  %                name        default   validation
  addParameter(o , 'l2'      , 0    ,    @isscalar);
  addParameter(o , 'maxiter' , 1000 ,    @isscalar);
  addParameter(o , 'tol'     , 1e-8 ,    @isscalar);
  parse(o, optcell{:});
  
  maxiter = o.Results.maxiter;
  l2      = o.Results.l2;
  tol     = o.Results.tol;

  num_tasks  = length(X);
  dimension = size(X{1}, 1);

  % initialize (can provide your own initialization)
  W0 = zeros(dimension, num_tasks);

  grad_flag = 0;

  Wz = W0;
  Wz_old = W0;

  t = 1;
  t_old = 0;

  iter = 1;
  gamma = 1;
  gamma_inc = 2;
  obj = zeros(maxiter,1);

  while iter < maxiter
    zeta = (t_old - 1) /t;

    Ws = (1 + zeta) * Wz - zeta * Wz_old;

    % compute function value and gradients of the search point
    [grad, Fs ]  = gradVal_eval(Ws);

    % the Armijo Goldstein line search
    while true
      if alpha>0
        Wzp = soslasso_projection(Ws - grad/gamma,alpha/gamma,lambda,group_arr,groups);
      else
        Wzp = Ws - grad/gamma;
      end
      Fzp = funVal_eval(Wzp);

      delta_Wzp = Wzp - Ws;
      nrm_delta_Wzp = norm(delta_Wzp, 'fro')^2;
      r_sum = (nrm_delta_Wzp);


      Fzp_gamma = Fs + sum(sum(delta_Wzp.* grad)) + gamma/2 * nrm_delta_Wzp;

      if (r_sum <=1e-20)
        grad_flag=1; % this shows that, the gradient step makes little improvement
        break;
      end

      if (Fzp <= Fzp_gamma)
        break;
      else
        gamma = gamma * gamma_inc;
      end
    end

    Wz_old = Wz;
    Wz = Wzp;


    if (grad_flag)
      break;
    end
    if alpha>0
      obj(iter) = Fzp + sos_eval(Wz,group_arr ,alpha,lambda);
    else
      obj(iter) = Fzp;
    end

    % convergence check.
    if iter>=2
      if (abs( obj(end) - obj(end-1) ) <= tol*obj(end-1))
        break;
      end
    end
    iter = iter + 1;
    t_old = t;
    t = 0.5 * (1 + (1+ 4 * t^2)^0.5);

  end

  W = Wzp;
  obj(iter+1:end) = [];
  % private functions

	function Wshr = soslasso_projection(W,lam,lambda,group_arr,groups)
		% step 1: perform soft thresholding
		X_soft = sign(W).*max(abs(W) - lambda*alpha,0);

		%step 2: perform group soft thresholding

		X_soft = [X_soft; zeros(1,size(X_soft,2))]; % for the dummy
		Xtemp = sum(X_soft.^2,2); %xtemp is now a vector
		Xtemp = sum(Xtemp(group_arr),2);
		Xtemp = sqrt(Xtemp);
		Xtemp = max(Xtemp - lam,0); % this is the multiplying factor
		Xtemp = Xtemp./(Xtemp + lam);
		Xtemp = Xtemp(groups);
		Xtemp = repmat(Xtemp,1,size(X_soft,2));
		Wshr = X_soft(1:end-1,:).*Xtemp;

	end


	function [grad_W, funcVal] = gradVal_eval(W)
		grad_W = zeros(dimension, num_tasks);
		lossValVect = zeros (1 , num_tasks);

		for ii = 1:num_tasks
			[ grad_W(:, ii), lossValVect(:, ii)] = unit_grad_eval( W(:, ii), X{ii}, Y{ii});
		end

		grad_W = grad_W + l2 * 2 * W;

		funcVal = sum(lossValVect) + l2*norm(W,'fro')^2;
	end

	function [funcVal] = funVal_eval (W)
		funcVal = 0;

			for ii = 1: num_tasks
				funcVal = funcVal + unit_funcVal_eval( W(:, ii), X{ii}, Y{ii});
			end

		funcVal = funcVal + l2 * norm(W,'fro')^2;
	end


% SOS regularizer value

	function [regval] = sos_eval(W,group_arr ,alpha,lambda)
		regval = 0;
		[n,p] = size(group_arr);
		Wtemp = [W;zeros(1,num_tasks)];

		for ii = 1 : n
			w = Wtemp(unique(group_arr(ii,:)), :);
			%w = w.^2;
			%regval = regval + alpha * sqrt(sum(w(:)));
      regval = regval + (alpha * norm(w,2));
		end
		regval = regval + alpha*lambda*norm(W(:),1);

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
% ------------------------------------------------------------------------
function opt = getDefaultField(data,field,default)
% ------------------------------------------------------------------------
  if isfield(data,field)
    opt = data.(field);
  else
    opt = default;
  end
end
