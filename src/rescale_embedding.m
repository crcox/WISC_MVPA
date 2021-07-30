% RESCALE_EMBEDDING Apply singular values to scale embedding dimensions.
%
%Embedding is achieved via singulat value decomposition. The orthonormal
%singular vectors can be rescaled to reflect their importance to the
%reconstruction by multiplying by the singular values.
%
% Syntax:  C = rescale_embedding(U, z)
%
% Inputs:
%    U - Matrix of orthonormal singular vectors.
%    z - Singular values corresponding to the vectors in C.
%
% Outputs:
%    C - The r most important orthonormal singular vectors.
%
% Example:
%    S = corr(randn(3, 100));
%    r = embedding_rank_at_criterion(S, 0.2);
%    [C0, z] = embed_similarity_matrix(S, r);
%    C = rescale_embedding(C0, z);
%
% Other m-files required: none
% Subfunctions: objective_function
% MAT-files required: none
%
% See also: EMBEDDING_RANK_AT_CRITERION RESCALE_EMBEDDING
%
% Author: Christopher R. Cox
% Department of Psychology, Louisiana State University
% email: chriscox@lsu.edu
% Website: https://faculty.lsu.edu/chriscox
% Last revision: 29-July-2021
function C = rescale_embedding(U, z)
    C = U * diag(sqrt(z));
end
