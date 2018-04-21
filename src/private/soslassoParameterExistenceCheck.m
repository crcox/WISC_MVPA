function [b, c] = soslassoParameterExistenceCheck(p)
% SOSLASSOPARAMETEREXISTENCECHECK Checks if valid combination of parameters
% was passed.
%
% SOS Lasso is a complicated optimization, and the interface to it
% continues to evolve. For example, sometimes it is more appropriate to
% consider the regularization parameters in terms of alpha and lambda,
% where alpha is a ratio controlling the relative importance of the
% grouping penalty vs. overall sparsity. The idea is to scale between group
% lasso and lasso. However, it is not well behaved when alpha is very small
% (approaching pure lasso). In order to actually implement lasso, absolute
% control over the grouping penalty (lamSOS) is required. The former
% parameterization scheme is refered to as "ratio parameters" and the
% latter is refered to as "absolute parameters".
%
% Critial to the SOS Lasso optimization is the pre-assignment of
% features into groups. This can be done either by specifying parameters
% that will generate groups "on the fly" in terms of diameter, amount of
% overlap, and shape. Alternatively, groups might be specified by some
% other means and stored in a file to be referenced. The former
% parameterization scheme is refered to as "parametric groups", and the
% latter is refered to as "premade groups".
%
% This function assesses which reguarization parameter scheme (ratio or
% absolute) and which group-defining scheme (parametric or premade) is
% being used.
%
% INPUT
%  p : the parameter structure
%
% OUTPUT
%  b : TRUE if a valid parameter scheme is detected, FALSE if not.
%  c : A structure identifying which scheme is present, one of either:
%       1.  'ratio parameters and premade groups'
%       2.  'absolute parameters and premade groups'
%       3.  'ratio parameters and parametric groups'
%       4.  'absolute parameters and parametric groups'
%  
    B = true(2,2);
    condition_combos = {
        'ratio parameters and premade groups'
        'absolute parameters and premade groups'
        'ratio parameters and parametric groups'
        'absolute parameters and parametric groups'
    };
    % Which Regularization parameters?
    % --------------------------------
    % 1. Using lamSOS and lamL1?
    B(1,1) = B(1,1) && isfield(p,'lamSOS') && ~isempty(p.lamSOS);
    B(1,1) = B(1,1) && isfield(p,'lamL1') && ~isempty(p.lamL1);
    % 2. Using alpha and lambda?
    B(2,1) = B(2,1) && isfield(p,'alpha') && ~isempty(p.alpha);
    B(2,1) = B(2,1) && isfield(p,'lambda') && ~isempty(p.lambda);
    
    % Which Group-defining parameters?
    % --------------------------------
    % 1. Using diameter, shape, and overlap?
    B(1,2) = B(1,2) && isfield(p,'diameter') && ~isempty(p.diameter) && isnumeric(p.diameter);
    B(1,2) = B(1,2) && isfield(p,'shape') && ~isempty(p.shape) && (ischar(p.shape) || iscellstr(p.shape));
    B(1,2) = B(1,2) && isfield(p,'overlap') && ~isempty(p.overlap);
    % 2. Using sosgroups?
    B(2,2) = B(2,2) && isfield(p,'sosgroups') && ~isempty(p.sosgroups);
    
    b = all(xor(B(1,:),B(2,:)));
    if b
        k = combvec([0,1],[0,1]);
        z = all(bsxfun(@eq,k,B(1,:)'));
        c = struct('id',find(z),'description',condition_combos{z});
    else
        c = struct('id',0,'description','INVALID');
    end
end
