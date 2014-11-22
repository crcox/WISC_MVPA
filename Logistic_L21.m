%% FUNCTION Logistic_L21
%  L21 Joint Feature Learning with Logistic Loss.
%
%% OBJECTIVE
%   argmin_{W,C} { sum_i^t (- sum(log (1./ (1+ exp(-X{i}*W(:, i) - Y{i} .* C(i)))))/length(Y{i}))
%			+ opts.rho_L2 * \|W\|_2^2 + rho1 * \|W\|_{2,1} }
%
%% INPUT
%   X: {n * d} * t - input matrix
%   Y: {n * 1} * t - output matrix
%   rho_L2: Regularizer on the loss function---does not need to be tuned; must be >0.
%     N.B. We use a loss function that behaves like logistic regression,
%     but is not exactly logistic regression. It depends on being
%     regularized to prevent driving values to infinity. We are maximizing
%
%% OUTPUT
%   W: model: d * t
%   funcVal: function value vector.

%
%% RELATED FUNCTIONS
%  Least_L21, init_opts

%% Code starts here
function [W, funcVal,iter] = Logistic_L21(X, Y, GroupInfo, rho_L2,lambda,alpha,opts)
	group_arr = GroupInfo.group_arr;
	groups = GroupInfo.groups;
	RepIndex = [1, GroupInfo.RepIndex+1]; % to account for bias term

	task_num  = length(X);
	for i = 1:task_num
		X{i} = X{i}';
	end

	% initialize options.
	dimension = length(RepIndex);
	funcVal = [];

	if isfield(opts,'init')
		if opts.init==1;
			Wz = [opts.C0;opts.W0];
		elseif opts.init==2;
			Wz = zeros(dimension, task_num); % initialize
		end
	end

	% this flag tests whether the gradient step only changes a little
	bFlag=0;

	Wz_old = Wz;

	t = 1;
	t_old = 0;
	iter = 0;
	gamma = 1;
	gamma_inc = 2;
	maxIter = opts.max_iter;
	tol = opts.tol;
	tFlag = opts.tFlag;
	while iter < maxIter
		zeta = (t_old - 1) /t;

		Ws = (1 + zeta) * Wz - zeta * Wz_old;

		% compute function value and gradients of the search point
		% This nested function accesses Ws, X, Y, and RepIndex.
		[gWs, Fs ]  = gradVal_eval(Ws,rho_L2);

		% the Armijo Goldstein line search scheme
		while true
			% This is the function that applies the SOS LASSO shrinkage.
			% Don't shrink the intercept term.
			x = Ws(2:end,:) - gWs(2:end,:)/gamma;
			Wzp = [Ws(1,:); soslasso_shrink_logistic(x,group_arr,groups,lambda/gamma,alpha)];

			% This nested function references X, Y, and RepIndex.
			Fzp = funVal_eval_PV(Wzp);
			% disp([sum(Wzp_1(:)),sum(Fzp(:))])

			delta_Wzp = Wzp - Ws;
			nrm_delta_Wzp = norm(delta_Wzp, 'fro')^2;
			r_sum = (nrm_delta_Wzp)/2;

			Fzp_gamma = Fs + sum(sum(delta_Wzp.* gWs)) + gamma/2 * nrm_delta_Wzp;

			if (r_sum <=1e-20)
			% This indicates that the gradient step results in virtually no
			% change.
				bFlag=1;
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

		funcVal = cat(1, funcVal, Fzp + nonsmooth_eval(Wz(2:end,:), lambda,alpha,group_arr));

		if (bFlag)
			break;
		end

		  % test stop condition.
		switch(tFlag)
		case 0
			if iter>=2
				if (abs( funcVal(end) - funcVal(end-1) ) <= tol)
					break;
				end
			end
		case 1 % This is the default.
			if iter>=2
				if (abs( funcVal(end) - funcVal(end-1) ) <= tol* funcVal(end-1))
					break;
				end
			end
		case 2
			if ( funcVal(end)<= tol)
				break;
			end
		  case 3
			if iter>=maxIter
				break;
			end
		end
		iter = iter + 1;
		t_old = t;
		t = 0.5 * (1 + (1+ 4 * t^2)^0.5);
	end
	W = Wzp;

	%% Nested Functions
	% N.B. Nested functions can access variables in the parent scope
	% without copying---they can also modify those variables, and the
	% consequences will be felt in the parent scope as well. In this case,
	% they are used basically to make sure the X does not get duplicated.
	function [gWs, Fs ]  = gradVal_eval(Ws,lambda)
	% compute gradient of the loss function at point Ws
		gWs = zeros(dimension, task_num);
		Fs = 0;
		for ii = 1:task_num
			gWs(:,ii) = -X{ii}(RepIndex,:)*Y{ii} + lambda*Ws(:,ii);
			Fs = Fs + (X{ii}(RepIndex,:)*Y{ii})'*Ws(:,ii) + lambda/2*norm(Ws(:,ii))^2;
		end
		Fs = -Fs;
	end

	function Fzp = funVal_eval_PV(Wzp)
    Fzp = 0;
    for ii = 1:task_num
      Fzp = Fzp + Y{ii}'*X{ii}(RepIndex,:)'*Wzp(:,ii);
    end
  Fzp = -Fzp;
	end
end
