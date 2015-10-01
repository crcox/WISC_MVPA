% Use a shrinkage estimator of the covariance matrix from a sample
%
% In:
% - X - #data points x #variables
% - optional -
%   - 'target' - shrinkage target:
%     - 'diagonal' (default)
%
% Out:
% - U - the estimated covariance matrix
% - lambdaStar - the optimum value of the tradeoff parameter between sample covariance matrix and shrinkage target
%
% Notes:
% - implements one of the Ledoit/Wolff shrinkage estimators, as described in Schaffer/Strimmer 2005 
% - should be easy to add more shrinkage targets
%
% History
% 2009 Mar 13 - fpereira@princeton.edu - created from previous code
%

function [U,lambdaStar] = cov_shrinkage(varargin)

%% argument processing

if nargin < 1
  fprintf('syntax: cov_shrinkage <data matrix> [''target'',<''diagonal''|>]\n'); return;
else
  X = varargin{1}; [n,p] = size(X); % #points x #variables
  target = 'diagonal';
  algorithm = 'safe';
  
  idx = 2;
  while idx <= nargin
    argval = varargin{idx}; idx = idx + 1;
    
    switch argval
     case {'target'}
      target = varargin{idx}; idx = idx + 1;
     case {'algorithm'}
      algorithm = varargin{idx}; idx = idx + 1;
     otherwise
      fprintf('warning: unrecognized option %s, it will be ignored\n',argval);
    end
  end
  
  switch target
   case {'diagonal'}
   otherwise
    fprintf('error: shrinkage target %s is not supported\n',target);return;
  end

  switch algorithm
   case {'safe','fast'}
   otherwise
    fprintf('error: algorithm has to be safe|fast\n');return;
  end

end


%%% compute estimate via the parts specified in appendix A of paper

% subtract the mean from each variable
Xmean = mean(X,1);
X = X - repmat(Xmean,n,1);

% compute:
% - S    (sample covariance, p x p)
% - Svar (variance of sample covariance)

W = zeros(p,p,n);

if 1
  % slow, clear version

  % compute W and Wbar
  Wbar = zeros(p,p);
  for e = 1:n
    W(:,:,e) = X(e,:)'*X(e,:);
    Wbar = Wbar + W(:,:,e);
  end
  Wbar = Wbar / n;
  S    = n*Wbar/(n-1);
  Svar = zeros(p,p);
  for e = 1:n
    Svar = Svar + (W(:,:,e) - Wbar).^2;
  end
  Svar = n * Svar / ((n-1)^3);

  %% using S and Svar, compute the estimate
  
  S_diagonal       = diag(S);
  S_nondiagonal    = S(find(~eye(p)));  
  Svar_nondiagonal = Svar(find(~eye(p)));  
  
  lambdaStar = sum(Svar_nondiagonal)/sum(S_nondiagonal.^2);

  U = lambdaStar * diag(S_diagonal) + (1-lambdaStar) * S;
else
end
  



