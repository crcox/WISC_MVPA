function [results,info] = learn_category_encoding(Y, X, Gtype, varargin)
  p = inputParser();
  addRequired(p  , 'Y'                         );
  addRequired(p  , 'X'                         );
  addRequired(p  , 'Gtype'                     );
  addParameter(p , 'lambda'         , []       );
  addParameter(p , 'alpha'          , 0        );
  addParameter(p , 'groups'         , {}       );
  addParameter(p , 'cvind'          , []       );
  addParameter(p , 'cvholdout'      , []       );
  addParameter(p , 'normalize'      , []       );
  addParameter(p , 'DEBUG'          , false    );
  addParameter(p , 'AdlasOpts'      , struct() );
  addParameter(p , 'SmallFootprint' , false    );
  parse(p, Y, X, Gtype, varargin{:});

  Y         = p.Results.Y;
  X         = p.Results.X;
  G         = p.Results.groups;
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
  
  % For cases like lasso, the group structure is irrelevant, but in the
  % spirit of running every analysis through the same function, noG defines
  % a single group for each subject that contains all voxels.
  noG = cell(1,numel(X));
  for i = 1:numel(noG)
    noG{i} = 1:size(X{i},2);
  end

  if isempty(LAMBDA)
    nlam = 1;
  else
    nlam = length(LAMBDA);
  end

  if isempty(ALPHA)
    nalpha = 1;
  else
    nalpha = length(ALPHA);
  end

  ncv = max(cvind{1});
  if isempty(holdout)
    cvset = 1:ncv;
  else
    cvset = holdout;
  end

  zz = zeros(ncv,nalpha,nlam);
  h1 = cell(size(X)); [h1{:}] = deal(zz);
  h2 = cell(size(X)); [h2{:}] = deal(zz);
  f1 = cell(size(X)); [f1{:}] = deal(zz);
  f2 = cell(size(X)); [f2{:}] = deal(zz);
  m1 = cell(size(X)); [m1{:}] = deal(zz);
  m2 = cell(size(X)); [m2{:}] = deal(zz);
  c1 = cell(size(X)); [c1{:}] = deal(zz);
  c2 = cell(size(X)); [c2{:}] = deal(zz);
  err1 = cell(size(X)); [err1{:}] = deal(zz);
  err2 = cell(size(X)); [err2{:}] = deal(zz);
  clear zz;
  
  cc = cell(ncv,nalpha,nlam);
  BzAll   = cell(size(X)); [BzAll{:}] = deal(cc);
  YzAll   = cell(size(X)); [YzAll{:}] = deal(cc);
  clear cc;

  % Set number of tasks (i.e., subjects)
  if iscell(X)
    t = size(X,1);
  else
    t = 1;
  end
  
  fprintf('%7s%5s %7s %7s %11s %11s %11s %11s %11s\n', '','subj','alpha','lambda','test err','train err','test diff','train diff','n vox')
	for i = cvset
    X = Xorig;
    
    test_set = cell(size(cvind));
    train_set = cell(size(cvind));
    for ii = 1:numel(cvind)
      test_set{ii}  = cvind{ii}==i;
      train_set{ii} = ~test_set{ii};
    end
    

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
      if isempty(ALPHA)
        alpha = nan(1);
      elseif length(ALPHA) > 1
        alpha = ALPHA(j);
      else
        alpha = ALPHA;
      end

      for k = 1:nlam
        if isempty(LAMBDA)
          lambda = nan(1);
        elseif length(LAMBDA) > 1
          lambda = LAMBDA(k);
        else
          lambda = LAMBDA;
        end

        if DEBUG
          Bz = cell(size(X));
          for ii = 1:numel(X);
            Bz{ii} = randn(size(X{ii},2),1);
          end
          info.message = 'DEBUG';
          info.iter    = 0;
        else
          switch Gtype
          case 'lasso'
            [Bz, info] = SOS_logistic( ...
              subsetAll(X, train_set), ...
              subsetAll(Y, train_set), ...
                     0, lambda(k),noG, ...
              'l2'      ,    0, ...
              'maxiter' , 1000, ...
              'tol'     , 1e-8, ...
              'W0'      ,   []);

          case 'soslasso'
            [Bz, info] = SOS_logistic( ...
              subsetAll(X, train_set), ...
              subsetAll(Y, train_set), ...
              alpha,       lambda,  G, ...
              'l2'      ,    0, ...
              'maxiter' , 1000, ...
              'tol'     , 1e-8, ...
              'W0'      ,   []);
          end
        end
        for ii = 1:numel(X);
          BzAll{ii}{i,j,k} = Bz{ii};
          YzAll{ii}{i,j,k} = X{ii}*Bz{ii};
          Yz = YzAll{ii}{i,j,k};
          h1{ii}(i,j,k)   = nnz( Y{ii}(test_set{ii})  & (Yz(test_set{ii})>0)  ) / nnz( Y{ii}(test_set{ii}));
          h2{ii}(i,j,k)   = nnz( Y{ii}(train_set{ii}) & (Yz(train_set{ii})>0) ) / nnz( Y{ii}(train_set{ii}));
          f1{ii}(i,j,k)   = nnz(~Y{ii}(test_set{ii})  & (Yz(test_set{ii})>0)  ) / nnz(~Y{ii}(test_set{ii}));
          f2{ii}(i,j,k)   = nnz(~Y{ii}(train_set{ii}) & (Yz(train_set{ii})>0) ) / nnz(~Y{ii}(train_set{ii}));
          err1{ii}(i,j,k) = nnz( Y{ii}(test_set{ii})  ~= (Yz(test_set{ii})>0 )) / length(Y{ii}(test_set{ii}));
          err2{ii}(i,j,k) = nnz( Y{ii}(train_set{ii}) ~= (Yz(train_set{ii})>0)) / length(Y{ii}(train_set{ii}));
        end
        k1 = sum(cellfun(@nnz,Bz))/t;

        if isempty(ALPHA)
          alpha_j = nan;
        else
          alpha_j = ALPHA(j);
        end
        
        if isempty(LAMBDA)
          lambda_k = nan;
        else
          lambda_k = LAMBDA(k);
        end

        for ii = 1:numel(X)
          fprintf('cv %3d: %3d |%6.2f |%6.2f |%10.2f |%10.2f |%10.2f |%10.2f |%10.2f |\n', ...
            i,ii,alpha_j,lambda_k,err1{ii}(i,j,k),err2{ii}(i,j,k),h1{ii}(i,j,k)-f1{ii}(i,j,k),h2{ii}(i,j,k)-f2{ii}(i,j,k),nnz(BzAll{ii}{i,j,k}));
        end

%        fprintf('Exit status -- %s (%d iterations)\n', info.message, info.iter);
      end % lamba loop
    end % alpha loop
	end % cv loop
  for ii = 1:numel(X)
    if ~SMALL
      results(ii).Bz = BzAll{ii}(cvset,:,:);
      results(ii).Yz = YzAll{ii}(cvset,:,:);
    end
    results(ii).h1   = h1{ii}(cvset,:,:);
    results(ii).h2   = h2{ii}(cvset,:,:);
    results(ii).f1   = f1{ii}(cvset,:,:);
    results(ii).f2   = f2{ii}(cvset,:,:);
    results(ii).m1   = m1{ii}(cvset,:,:);
    results(ii).m2   = m2{ii}(cvset,:,:);
    results(ii).c1   = c1{ii}(cvset,:,:);
    results(ii).c2   = c2{ii}(cvset,:,:);
    results(ii).err1 = err1{ii}(cvset,:,:);
    results(ii).err2 = err2{ii}(cvset,:,:);
    results(ii).subject = ii;
    results(ii).note = 'Matrix dims = 1:cv 2:alpha 3:lambda';
  end
end

function X = doNormalization(X, train_set) %#ok<DEFNU>
  mm = ones(n,1)*mean(X(train_set,:),1);
  ss = ones(n,1)*std(X(train_set,:),1);
  if all(X(:,end)==1)
    ss(:,d) = ones(n,1);
  end
  X = (X-mm)./ss;
end
