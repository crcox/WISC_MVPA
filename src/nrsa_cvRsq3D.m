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
  ZSCORE = false;
  if ZSCORE
    x = zscore(x);
  end
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

    % If there are any NaN values at this stage, that indicates an
    % ideosyncratic filter applies to one of the conditions. All typical
    % conditions can be handled in one shot, but ideosyncratic conditions will
    % be handled individually.
    z = all(squeeze(any(isnan(pp),1)),2);
    for i_r = 1:r
      if ZSCORE
        pp(:, ~z, i_r) = zscore(pp(:, ~z, i_r));
      end
      y(~z, i_r) = 1 - (mean(bsxfun(@minus, pp(:, ~z, i_r), x(:, i_r)).^2) ./ xvar(i_r));
    end
    if any(z)
      ix = find(z);
      for i_nan = 1:nnz(z)
        for i_r = 1:r
          z_nan = isnan(pp(:, ix(i_nan), i_r));
          if ZSCORE
            pp(~z_nan, ix(i_nan), i_r) = zscore(pp(~z_nan, ix(i_nan), i_r));
          end
          y(ix(i_nan), i_r) = 1 - (mean(bsxfun(@minus, pp(~z_nan, ix(i_nan), i_r), x(~z_nan, i_r)).^2) ./ xvar(i_r));
        end
      end
    end

  elseif ndims(p) == 2
    if ZSCORE
      p = zscore(p);
    end
    y = 1 - (mean((p - x).^2) ./ xvar);
  end
end
