function [Xhat,C,iter] = overlap_2stage_CMU(Y,Xo,X,G,group_arr,groups, lambda,gamma)

% function to perform the overlapping SGL optimization, with both LS and LOGIT loss
% INPUTS
% loss    = type of loss function. 0 = least squares, 1 = logistic
% Y       = T X 1 cell array of observations . T  =number of tasks
% X       = T X 1 cell array of data
% Xo      = T X 1 cell array of replicated data matrices
% G       = cell array of groups
% group_arr = output from replication step
% lambda  = regularizer
% OUTPUTS
% Xhat    = debiased output
% C       = the bias term in the output
% 
% CODE REQUIRES MALSAR PACKAGE
%
% Nikhil Rao
% 3/17/13

debias = false;

[W, C, ~,iter] = Logistic_L21_CMU(Xo, Y, lambda, gamma, group_arr, groups);

% W is the output matrix
%we now need to combine overlapping groups
n = zeros(length(G),1);
for iii = 1:length(G)
    temp = G{iii};
    n(iii) = max(temp);
end
n = max(n);
T = length(Y);
Xhat = zeros(n,T);

% identify whether a dummy variable exists and chuck it
dummy = max(max(group_arr));
mask = (group_arr == dummy);
isdummy = 0;
if sum(sum(mask))>1
    isdummy = 1;
end
for ii = 1:length(G)
    t = G{ii};
    s = group_arr(ii,:);
    if isdummy == 1
       dummyind = find(s == dummy);
       s(dummyind) = [];
    end
    Xhat(t,:) = Xhat(t,:) + W(s,:);
end

% debias the solution using the X and Y and Xhat
if debias
    inds = cell(T,1);
    temp = Xhat;
    Xhat = zeros(n,T);
    for ii = 1:T
        idx = find(temp(:,ii)~=0);
        if ~isempty(idx);
            inds{ii} = idx;
            Xtemp = X{ii};
            Xtemp = Xtemp(:,idx);
            Bhat = Xtemp\Y{1};
            Xhat(idx,ii) = Bhat;
        end
    end
end
   
end




