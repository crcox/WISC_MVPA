% EMBED_RANK_AT_CRITERION Determine rank to achieve reconstruction
% accuracy.
%
%Following a singular value decomposition of a matrix S, S can be recomposed
%from the singular vectors and singular values. We seek a low-rank embedding of
%S that includes only the r most important singular vectors and values. As r
%increases, the fidelity of the reconstruction will increase.
%
%This function will determine the smallest r that achieves a specified
%criterion tau. Tau is specified in units of the objective function (which
%monotonically decreases as r increases).
%
% Syntax:  [r] = embed_rank_at_criterion(S, tau)
%
% Inputs:
%    S   - A symmetric matrix of pairwise similarity values.
%    tau - A target value for the loss function.
%
% Outputs:
%    r - The rank of the smallest embedding that achieves the criterion defined
%        by tau.
%
% Example:
%    S = corr(randn(3, 100));
%    r = embed_rank_at_criterion(S, 0.2)
%    C = embed_similarity_matrix(S, r);
%
% Other m-files required: none
% Subfunctions: objective_function
% MAT-files required: none
%
% See also: EMBED_SIMILARITY_MATRIX RESCALE_EMBEDDING
%
% Author: Christopher R. Cox
% Department of Psychology, Louisiana State University
% email: chriscox@lsu.edu
% Website: https://faculty.lsu.edu/chriscox
% Last revision: 29-July-2021
function r = embed_rank_at_criterion(S, tau)
    [U, Z, ~] = svd(S);
    z = diag(Z);
    for r = 1:size(U, 2);
        C = U(:, 1:r) * diag(sqrt(z(1:r)));
        if obj_func(S, C * C') <= tau
            break
        end
    end
end

function loss = objective_function(y, p)
    loss = norm(y - p, 'fro') / norm(y, 'fro');
end
