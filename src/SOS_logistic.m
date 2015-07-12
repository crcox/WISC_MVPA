function [W, obj] = SOS_logistic(X, Y,alpha,lambda,G,varargin)

% function solves
% min \sum_t log(1 + exp(-Yt.*XtW(:,t))) + lambda Omega(W)
% 
%
% INPUT
% X         : {n x d} length t cell array
% Y         : {n x 1} length t output
% alpha     : Ratio of L1 to SOS Lasso penalty importance.
%             0 --> lasso
%             1 --> group lasso
% lambda    : Scales overall severity of regularization.
%             SOS Lasso penalty = lambda * alpha
%             L1 penalty = lambda * (1 - alpha)
% G         : Cell array of group assignments for each subject. Each row
%             corresponds to a a group, and each column corresponds to a 
%             subject. Values in each cell are indexes that refer to 
%             columns in X.
% -- optional --
% l2        : an additional l2 regularizer for smoothing the solution.
%             Default is 0.
% maxiter   : max iterations. Default is 1000.
% tol       : Tolerance for determining convergence. Default is 1e-8.
%
% OUTPUT
% W: model: d x t

%%%%%%%%%
% parts of this code have been adapted from the MALSAR package
% http://www.public.asu.edu/~jye02/Software/MALSAR/
%%%%%%%%%

  p = inputParser;
  %% Parse function inputs
  %                name         default     validation
  addRequired(p  , 'X'                                );
  addRequired(p  , 'Y'                                );
  addRequired(p  , 'alpha'                , @isscalar );
  addRequired(p  , 'lambda'               , @isscalar );
  addRequired(p  , 'G'                    , @iscell   );
  addParameter(p , 'l2'      , 0          , @isscalar );
  addParameter(p , 'maxiter' , 1000       , @isscalar );
  addParameter(p , 'tol'     , 1e-8       , @isscalar );
  addParameter(p , 'W0'      , []         , @iscell   );
  parse(p, X, Y, alpha, lambda, G, varargin{:});

  X         = p.Results.X;
  Y         = p.Results.Y;
  alpha     = p.Results.alpha;
  lambda    = p.Results.lambda;
  G         = p.Results.G;
  maxiter   = p.Results.maxiter;
  l2        = p.Results.l2;
  tol       = p.Results.tol;
  W0        = p.Results.W0;
  
  % Setup group indexes
  [Gc, ix]  = commonGrouping(G);
  group_arr = group2mat(Gc);
  groups    = group2lab(Gc);
  for j = 1:length(X);
    X{j} = [X{j},zeros(size(X{j},1),1)];
    X{j} = X{j}(:,ix(:,j));
    X{j} = X{j}';
    if ~isempty(W0)
      W0{j} = [W0{j}(:); 0];
      W0{j} = W0{j}(ix(:,j));
    end
  end
  
  % Set lamL1 and lamSOS
  [lamSOS, lamL1] = ratio2independent(alpha, lambda);
  % Equivalent to:
  %   lamSOS = lambda * alpha;
  %   lamL1  = lambda * (1 - alpha);

  num_tasks = numel(X);
  dimension = length(groups);

  % initialize (can provide your own initialization)
  if isempty(W0)
    W0 = zeros(dimension, num_tasks);
  else
    W0 = cell2mat(W0);
  end

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
      if lamSOS>0
        Wzp = soslasso_projection(Ws - grad/gamma,lamSOS/gamma,lamL1,group_arr,groups);
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
    if lamSOS>0
      obj(iter) = Fzp + sos_eval(Wz,group_arr ,lamSOS,lamL1);
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
      
  W = combineOverlappingWeights(W,G);
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
    Xtemp = repmat(Xtemp,1,size(X,2));
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
