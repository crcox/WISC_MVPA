% EMBED_SIMILARITY_MATRIX Return embedding of specified rank
%
%Embedding is achieved via singulat value decomposition. The returned U matrix
%contain the orthonormal singular vectors, and the vector z contains the
%singular values which can be used to scale each vector relative to the amount
%of variance it accounts for in S.
%
% Syntax:  [C, z] = embed_similarity_matrix(S, r)
%
% Inputs:
%    S - A symmetric matrix of pairwise similarity values.
%    r - The desired rank of the embedding.
%
% Outputs:
%    C - The r most important orthonormal singular vectors.
%    z - Singular values corresponding to the vectors in C.
%
% Example:
%    S = corr(randn(3, 100));
%    r = embedding_rank_at_criterion(S, 0.2);
%    [U, z] = embed_similarity_matrix(S, r);
%    C = rescale_embedding(U, z);
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
function [U, z] = embed_similarity_matrix(S, r)
  [U, Z, ~] = svd(S);
  z = diag(Z);
  U = U(:, 1:r);
  z = z(1:r);
end
