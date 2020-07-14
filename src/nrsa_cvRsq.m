function y = nrsa_cvRsq3D(x, p, xvar)
% NRSA_CVRSQ Cross-validated R^2 (Q^2)
%
% INPUT
% x    : the target matrix
% p    : the predicted matrix
% xvar : (optional) fixed the denominator
%
% OUTPUT
% y : error term
  if nargin < 3
    xvar = var(x);
  end
  [m, r] = size(x);
  if ndims(p) == 3
    m_dim = find(size(p) == m);
    r_dim = find(size(p) == r);
    n_dim = find(~ismember(1:3, [m_dim, r_dim]));
    pp = permute(p, [m_dim, n_dim, r_dim]);
    n = size(p, n_dim);
    y = zeros(n, r);
    for i_r = 1:r
      y(:, i_r) = 1 - (mean(bsxfun(@minus, p, x(:, i_r)).^2) ./ xvar);
    end

  elseif ndims(p) == 2
    y = 1 - (mean((p - x).^2) ./ xvar);
  end
end
