function [Xhat,Wndb,W,iter] = overlap_2stage(Y,X,GroupInfo,lambda, alpha, opts)
% function to perform the Sparse Overlapping Sets Lasso (SOS Lasso)
% optimization, with either least squares or LOGIT loss.
%
% INPUTS
% loss    = type of loss function. 0 = least squares, 1 = logistic
% Y       = T X 1 cell array of observations . T  =number of tasks
% X       = T X 1 cell array of data
% G       = cell array of groups
% group_arr = output from replication step
% lambda  = regularizer
%
% OUTPUTS
% Xhat    = debiased output
% C       = the bias term in the output
%
% CODE REQUIRES MALSAR PACKAGE
%
% Nikhil Rao and Chris Cox
% 4/12/2014
%
	G = GroupInfo.G;
	group_arr = GroupInfo.group_arr;

	[W, ~,iter] = Logistic_L21(X,Y,GroupInfo,10,lambda,alpha,opts);
	% W contains the solutions for each subject---there will be as many columns as
	% there are cells in X and Y. In the case of overlapping groups, there will be
	% as many rows as there are voxels *in the replicated space*.

	m = length(G);
	n = zeros(m,1);
	for k = 1:m
    n(k) = max(G{k});
	end
	n = max(n)+1; % to account for bias weight
	T = length(Y);
	Xhat = zeros(n,T);

	% Dummy variables are used to fill non-cortical space for each subject.
	% Identify whether a dummy variable exists and chuck it
  for k = 1:m
    t = G{k};
    s = group_arr(k,:)+1; % to account for bias weight
    s = s(~isnan(s));
    Xhat(t,:) = Xhat(t,:) + W(s,:);
	end
	Xhat(1,:) = W(1,:);

	% The solution before debiasing.
	Wndb = Xhat;

	%debias the solution using the X and Y and Xhat
	% Xnew = cell(T,1);
	inds = cell(T,1);
	temp = Xhat;
	Xhat = zeros(n,T);
	% totinds =[];
	warning('off','stats:glmfit:IllConditioned');
	warning('off','stats:glmfit:PerfectSeparation');
	warning('off','stats:glmfit:IterationLimit');
	for ii = 1:T
    z = temp(:,ii) ~= 0;
    z(1) = false;
    if any(z);
      inds{ii} = find(z);
      Xdeb = X{ii}(:,z);
      Bhat = glmfit(Xdeb, Y{ii}>0, 'binomial');
      z(1) = true;
      Xhat(z,ii) = Bhat;
    end
	end
	warning('on','stats:glmfit:IllConditioned');
	warning('on','stats:glmfit:PerfectSeparation');
	warning('on','stats:glmfit:IterationLimit');
end
