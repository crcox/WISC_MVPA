function y = nrsa_corr3D(x, p)
% NRSA_CORR Correlation-base performance metric for Network RSA
%
% Computes the Spearman correlation between a target similarity matrix and
% a predicted similarity matrix (lower triangles).
%
% INPUT
% x : the target matrix
% p : the predicted matrix
%
% OUTPUT
% y : error term (the loss value)
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
      y(~z, i_r) = (zscore(x(:, i_r), 1)' * zscore(pp(:, ~z, i_r), 1)) ./ m;
    end
    if any(z)
      ix = find(z);
      for i_nan = 1:nnz(z)
        for i_r = 1:r
          z_nan = isnan(pp(:, ix(i_nan), i_r));
          y(ix(i_nan), i_r) = (zscore(x(~z_nan, i_r), 1)' * zscore(pp(~z_nan, ix(i_nan), i_r), 1)) ./ m;
        end
      end
    end

  elseif ndims(p) == 2
    y = zeros(1, r);
    for i_r = 1:r
      y(:, i_r) = (zscore(x(:, i_r), 1)' * zscore(p(:, i_r), 1)) ./ m;
    end
  end
end
