
function [W,obj] = SOG_least_squares(X, Y,lambda,lam1,groups,group_arr,varargin)

% function solves
% min \sum_t || Yt - XtW(:,t) ||^2 + lambda Omega(W)
% 
%
% INPUT
% X         : {n x d} length t cell array
% Y         : {n x 1} length t output
% lambda    : regularizer for the SOSlasso penalty
% lam1      : regularizer for the L1 part of soslasso penalty
% groups    : group assignment (see related function)
% group_arr : group matrix (see related functions)
% optional:
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



% optional parameter defaults
X = multi_transpose(X);

l2 = 0;
maxiter = 1000;
tol = 1e-8;

% check for optional parameters entered by the user
if (rem(length(varargin),2)==1)
    error('Optional parameters should be in pairs');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'maxiter'
                maxiter = varargin{i+1};
            case 'tol'
                tol = varargin{i+1};
            case 'l2'
                l2 = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end



num_tasks  = length (X);
dimension = size(X{1}, 1);

XY = cell(num_tasks, 1);
W0_prep = [];
for t_idx = 1 : num_tasks
    XY{t_idx} = X{t_idx}*Y{t_idx};
    W0_prep = cat(2, W0_prep, XY{t_idx});
end


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
    
    alpha = (t_old - 1) /t;
    
    Ws = (1 + alpha) * Wz - alpha * Wz_old;
    
    % compute function value and gradients of the search point
    grad  = gradVal_eval(Ws);
    Fout   = funVal_eval (Ws);
    
    % line search
    while true
        
        Wzp = soslasso_projection(Ws - grad/gamma,lambda/gamma,lam1,group_arr,groups);
        
        Fzp = funVal_eval(Wzp);
        
        delta_Wzp = Wzp - Ws;
        r_sum = norm(delta_Wzp, 'fro')^2;

        Fzp_gamma = Fout + sum(sum(delta_Wzp.* grad))...
            + gamma/2 * norm(delta_Wzp, 'fro')^2;
        
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
    
    obj(iter) = Fzp + sos_eval(Wz,group_arr ,lambda,lam1);
    
    if (grad_flag)
        break;
    end
    
    % convergence check
    
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

    function Wshr = soslasso_projection(W,lam,lam1,group_arr,groups)
        
        
        % step 1: perform soft thresholding
        X_soft = sign(W).*max(abs(W) - lam1*lambda,0);
        
        %step 2: perform group soft thresholding
      
        X_soft = [X_soft; zeros(1,size(X_soft,2))]; % for the dummy
        Xtemp = sum(X_soft.^2,2); %xtemp is now a vector
        Xtemp = sum(Xtemp(group_arr),2);
        Xtemp = sqrt(Xtemp);
        Xtemp = max(Xtemp - lam,0); % this is the multiplying factor
        Xtemp = Xtemp./(Xtemp + lam);
        Xtemp = Xtemp(groups);
        Xtemp = repmat(Xtemp,1,size(X,2));
        Wshr = X_soft(1:end-1,:).*Xtemp;
        
    end



% gradient of the loss function + l2 regularizer (if present)
    function [grad_W] = gradVal_eval(W)
        
        grad_W = zeros(dimension, num_tasks);
        for i = 1:num_tasks
            grad_W(:,i) = X{i}*(X{i}' * W(:,i)-Y{i});
        end
        
        grad_W = grad_W + l2*2*W;
    end

% loss function value + l2 regularizer (if present)
    function [funcVal] = funVal_eval(W)
        funcVal = 0;
        
        for i = 1: num_tasks
            funcVal = funcVal + 0.5 * norm (Y{i} - X{i}' * W(:, i))^2;
        end
        
        funcVal = funcVal + l2 * norm(W,'fro')^2;
    end


% SOS regularizer value

    function [regval] = sos_eval(W,group_arr ,lambda,lam1)
        regval = 0;
        [n,p] = size(group_arr);
        Wtemp = [W;zeros(1,num_tasks)];
        
        for i = 1 : n
            w = Wtemp(unique(group_arr(i,:)), :);
            w = w.^2;
            regval = regval ...
                + lambda * sqrt(sum(w(:)));
        end
        regval = regval + lambda*lam1*norm(W(:),1);
        
    end
end