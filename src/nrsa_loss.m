function y = nrsa_loss(x,p)
% NRSA_LOSS Loss function for Network RSA
%
% INPUT
% x : the target matrix
% p : the predicted matrix
%
% OUTPUT
% y : error term (the loss value)
    y = norm(x - p,'fro')/norm(x,'fro');
end