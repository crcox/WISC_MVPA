function [results,info] = learn_category_encoding(Y, X, Gtype, varargin)
  p = inputParser();
  addRequired(p  , 'Y'                         );
  addRequired(p  , 'X'                         );
  addRequired(p  , 'Gtype'                     );
  addParameter(p , 'lambda'         , []       );
  addParameter(p , 'alpha'          , 0        );
  addParameter(p , 'cvind'          , []       );
  addParameter(p , 'cvholdout'      , []       );
  addParameter(p , 'normalize'      , []       );
  addParameter(p , 'DEBUG'          , false    );
  addParameter(p , 'AdlasOpts'      , struct() );
  addParameter(p , 'SmallFootprint' , false    );
  parse(p, S, X, Gtype, varargin{:});

  S         = p.Results.S;
  X         = p.Results.X;
  Gtype     = p.Results.Gtype;
  LAMBDA    = p.Results.lambda;
  ALPHA     = p.Results.alpha;
  cvind     = p.Results.cvind;
  holdout   = p.Results.cvholdout;
  normalize = p.Results.normalize;
  DEBUG     = p.Results.DEBUG;
  options   = p.Results.AdlasOpts;
  SMALL     = p.Results.SmallFootprint;

  Xorig = X;

  if isempty(lambda)
    nlam = 1;
  else
    nlam = length(lambda);
  end

  if isempty(alpha)
    nalpha = 1;
  else
    nalpha = length(alpha);
  end

  ncv = max(cvind);
  if isempty(holdout)
    cvset = 1:ncv;
  else
    cvset = holdout;
  end

  h1      = cell(ncv,nlam,nalpha);
  f1      = cell(ncv,nlam,nalpha);
  d1      = cell(ncv,nlam,nalpha);
  d2      = cell(ncv,nlam,nalpha);
  dp1     = cell(ncv,nlam,nalpha);
  dp2     = cell(ncv,nlam,nalpha);
  err1    = cell(ncv,nlam,nalpha);
  err2    = cell(ncv,nlam,nalpha);
  nz_rows = cell(ncv,nlam,nalpha);
  BzAll   = cell(ncv,nlam,nalpha);
  if nlam > 1 || nalpha > 1
    nz_rows = squeeze(mat2cell(nz_rows, ncv, d, ones(1,nlam), ones(1,nalpha)));
  end

  % Set number of tasks (i.e., subjects)
  if iscell(X)
    t = size(X,1);
  else
    t = 1;
  end

  fprintf('%8s%6s%11s %11s  %11s  %11s  %11s  %11s  %11s %11s  \n', '','lambda','alpha','test err','train err','test diff','cor train','n vox')
  for i = cvset
    X = Xorig;

    CVsize = nnz(cvind==i);
    test_set  = cvind==i;
    train_set = ~test_set;
    fprintf('cv %3d: ', i)

    for ii = 1:numel(X);
      switch normalize
        case 'zscore_train'
          mm = mean(X{ii}(train_set,:),1);
          ss = std(X{ii}(train_set,:),0,1);
        case 'zscore'
          mm = mean(X{ii},1);
          ss = std(X{ii},0,1);
        case 'stdev'
          mm = 0;
          ss = std(X{ii},0,1);
        case '2norm'
          mm = mean(X{ii},1);
          ss = norm(X{ii});
        otherwise
          mm = 0;
          ss = 1;
      end
      X{ii} = bsxfun(@minus,X{ii}, mm);
      X{ii} = bsxfun(@rdivide,X{ii}, ss);
    end

    for j = 1:nalpha
      if isempty(alpha)
        alpha = nan(1);
      end

      for k = 1:nlam
        if isempty(lambda)
          lambda = nan(1);
        end

        if DEBUG
          Bz = randn(d,t);
          info.message = 'DEBUG';
          info.iter    = 0;
        else
          switch Gtype
          case 'lasso'
            [Bz, info] = SOG_logistic(subsetAll(X, train_set), subsetAll(Y, train_set), ...
                                      0, lambda(k), [], [], options);

          case 'soslasso'

            [Bz, info] = SOG_logistic(subsetAll(X, train_set), subsetAll(Y, train_set), ...
                                      alpha(j), lambda(k), ...
                                      groups, group_arr, options);
          end
        end

        BzAll{i,j,k} = Bz;

        k1 = nnz(Bz)/t;
        Yz = X*Bz;

        SzAll{i,j,k} = Sz;

        %store results
        if nlam > 1 || nalpha > 1
          nz_rows{j,k}(i,:) = any(Bz,2);
        else
          nz_rows(i,:) = any(Bz,2);
        end

        % Comparison to true Y
        h1(i,j,k)   = nnz( Y(test_set)  & Yz(test_set)  / nnz( Y(test_set));
        h2(i,j,k)   = nnz( Y(train_set) & Yz(train_set) / nnz( Y(train_set));
        f1(i,j,k)   = nnz(~Y(test_set)  & Yz(test_set)  / nnz(~Y(test_set));
        f2(i,j,k)   = nnz(~Y(train_set) & Yz(train_set) / nnz(~Y(train_set));
        d1(i,j,k)   = h1(i,j,k) - f1(i,j,k);
        d2(i,j,k)   = h2(i,j,k) - f2(i,j,k);
        dp1(i,j,k)  = norminv(h1(i,j,k)) - norminv(f1(i,j,k));
        dp2(i,j,k)  = norminv(h2(i,j,k)) - norminv(f2(i,j,k));
        err1(i,j,k) = nnz( Y(test_set)  == Yz(test_set)  ) / length(Y);
        err2(i,j,k) = nnz( Y(train_set) == Yz(train_set) ) / length(Y);

        if isempty(lambda)
          lambda_j = nan;
        else
          lambda_j = lambda(j);
        end
        if isempty(alpha)
          alpha_k = nan;
        else
          alpha_k = alpha(k);
        end

        fprintf('%6.2f | %6.2f | %10.2f | %10.2f | %10.2f | %10.2f | %10.2f | %10.2f | %10d\n', ...
          lambda_j,alpha_k,err1(i,j,k),err2(i,j,k),p1(i,j,k),p2(i,j,k),cor1(i,j,k),cor2(i,j,k),k1);

        fprintf('Exit status -- %s (%d iterations)\n', info.message, info.iter);
      end % lamba loop
    end % alpha loop
  end % cv loop
  if ~SMALL
    results.Bz = UzAll;
    results.Sz = SzAll;
  end
  if ~iscell(nz_rows);
    results.nz_rows = nz_rows(cvset,:,:);
  else
    results.nz_rows = nz_rows;
  end
  results.h1   = t1(cvset,:,:);
  results.h2   = t2(cvset,:,:);
  results.f1   = f1(cvset,:,:);
  results.f2   = f2(cvset,:,:);
  results.d1   = d1(cvset,:,:);
  results.d2   = d2(cvset,:,:);
  results.dp1  = dp1(cvset,:,:);
  results.dp2  = dp2(cvset,:,:);
  results.err1 = err1(cvset,:,:);
  results.err2 = err2(cvset,:,:);
  results.iter = info.iter;
end % learn_similarity_encoding

function X = doNormalization(X, train_set)
  mm = ones(n,1)*mean(X(train_set,:),1);
  ss = ones(n,1)*std(X(train_set,:),1);
  if all(X(:,end)==1)
    ss(:,d) = ones(n,1);
  end
  X = (X-mm)./ss;
end
