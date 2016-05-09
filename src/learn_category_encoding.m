function [results,info] = learn_category_encoding(Y, X, regularization, varargin)
  p = inputParser();
  addRequired(p  , 'Y'                         );
  addRequired(p  , 'X'                         );
  addRequired(p  , 'regularization'                     );
  addParameter(p , 'lambda'         , []       );
  addParameter(p , 'alpha'          , 0        );
  addParameter(p , 'groups'         , {}       );
  addParameter(p , 'cvind'          , []       );
  addParameter(p , 'cvholdout'      , []       );
  addParameter(p , 'normalize'      , []       );
  addParameter(p , 'DEBUG'          , false    );
  addParameter(p , 'debias'         , false    );
  addParameter(p , 'AdlasOpts'      , struct() );
  addParameter(p , 'SmallFootprint' , false    );
  addParameter(p , 'PermutationTest' , false    );
  parse(p, Y, X, regularization, varargin{:});

  Y         = p.Results.Y;
  X         = p.Results.X;
  G         = p.Results.groups;
  regularization     = p.Results.regularization;
  LAMBDA    = p.Results.lambda;
  ALPHA     = p.Results.alpha;
  cvind     = p.Results.cvind;
  holdout   = p.Results.cvholdout;
  normalize = p.Results.normalize;
  DEBUG     = p.Results.DEBUG;
  DEBIAS    = p.Results.debias;
  options   = p.Results.AdlasOpts;
  SMALL     = p.Results.SmallFootprint;
  PermutationTest = p.Results.PermutationTest;

  Xorig = X;

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

  if isempty(holdout)
    cvset = 1:max(cvind{1});
  else
    cvset = holdout;
  end
  ncv = numel(cvset);

  % Set number of tasks (i.e., subjects)
  if iscell(X)
    t = numel(X);
  else
    t = 1;
  end

  % Define results struct
  results.Wz = [];
  results.Wix = [];
  results.Yz = [];
  results.Wnz = [];
  results.nvox = [];
  results.subject = [];
  results.cvholdout = [];
  results.finalholdout = [];
  results.alpha = [];
  results.lambda = [];
  results.diameter = [];
  results.overlap = [];
  results.shape = [];
  results.nt1  = [];
  results.nt2  = [];
  results.nd1  = [];
  results.nd2  = [];
  results.h1   = [];
  results.h2   = [];
  results.f1   = [];
  results.f2   = [];
  results.err1 = [];
  results.err2 = [];

  % Preallocate
  results(t*ncv*nlam*nalpha).Wz = [];

  % Permute if requested
  fprintf('PermutationTest: %d\n', PermutationTest);
  if PermutationTest
    n = size(Y,1);
    fprintf('Permuting %d rows of Y.\n', n);
    fprintf('First 10 rows of Y, before shuffling.\n')
    disp(Y(1:10,:))
    permix = randperm(n);
    if iscell(Y)
      for i = 1:numel(Y)
        y = Y{i};
        Y{i} = y(permix);
      end
    else
      Y = Y(permix);
    end
    fprintf('First 10 rows of Y, after shuffling.\n')
    disp(Y(1:10,:))
  end

  fprintf('%7s%5s %7s %7s %11s %11s %11s %11s %11s\n', '','subj','alpha','lambda','test err','train err','test diff','train diff','n vox')
  iii = 0;
  for i = 1:ncv
    icv = cvset(i);
    X = Xorig;

    test_set = cell(size(cvind));
    train_set = cell(size(cvind));
    for iSubject = 1:numel(cvind)
      test_set{iSubject}  = cvind{iSubject}==icv;
      train_set{iSubject} = ~test_set{iSubject};
    end

    for iSubject = 1:numel(X);
      switch lower(normalize)
        case 'zscore_train'
          mm = mean(X{iSubject}(train_set,:),1);
          ss = std(X{iSubject}(train_set,:),0,1);
        otherwise
          % All other cases already handled.
          mm = 0;
          ss = 1;
      end
      X{iSubject} = bsxfun(@minus,X{iSubject}, mm);
      X{iSubject} = bsxfun(@rdivide,X{iSubject}, ss);
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
          Wz = cell(size(X));
          for iSubject = 1:numel(X);
            Wz{iSubject} = randn(size(X{iSubject},2),1);
          end
          info.message = 'DEBUG';
          info.iter    = 0;
        else
          switch lower(regularization)
          case 'lasso_glmnet'
            [Wz, info] = lasso_glmnet( ...
              subsetAll(X, train_set), ...
              subsetAll(Y, train_set), ...
                       1, lambda,...
              'cvind'   , subsetAll(cvind, train_set),...
              'maxiter' , 1000, ...
              'tol'     , 1e-8, ...
              'W0'      ,   []);

          case 'lasso'
            [Wz, info] = SOS_logistic( ...
              subsetAll(X, train_set), ...
              subsetAll(Y, train_set), ...
                       1, lambda, G, ...
              'l2'      ,    0, ...
              'maxiter' , 1000, ...
              'tol'     , 1e-8, ...
              'W0'      ,   []);

          case 'soslasso'
            [Wz, info] = SOS_logistic( ...
              subsetAll(X, train_set), ...
              subsetAll(Y, train_set), ...
              alpha,       lambda,  G, ...
              'l2'      ,    0, ...
              'maxiter' , 1000, ...
              'tol'     , 1e-8, ...
              'W0'      ,   []);
          end
          if DEBIAS
            % Depends on glmnet
            Wz = ascell(Wz);
            opts_debias = glmnetSet(struct('alpha',0,'intr',0));
            for iSubject = 1:size(Wz)
              z = Wz{iSubject} ~= 0;
              if nnz(z) > 0
                cv = cvind{iSubject};
                cv = cv(cv~=icv);
                cv(cv>icv) = cv(cv>icv) - 1;

                debiasObj = cvglmnet(X{iSubject}(train_set{iSubject},z),Y{iSubject}(train_set{iSubject}), ...
                           'binomial',opts_debias,'class',max(cv),cv);
                opts_debias.lambda = debiasObj.lambda_min;
                debiasObj = glmnet(X{iSubject}(train_set{iSubject},z),Y{iSubject}(train_set{iSubject}),'binomial',opts_debias);
                Wz{iSubject}(z) = debiasObj.beta;
              end
            end
          end
        end

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

        for iSubject = 1:numel(X);
          iii = iii + 1;

          wz = Wz{iSubject};
          ix = find(wz);
          nv = numel(wz);
          wnz = nnz(wz);
          wz = wz(ix);

          y = Y{iSubject};
          yz = X{iSubject} * Wz{iSubject};
          h1   = nnz( y(test_set{iSubject})>0  & (yz(test_set{iSubject})>0)  );
          h2   = nnz( y(train_set{iSubject})>0 & (yz(train_set{iSubject})>0) );
          f1   = nnz(~y(test_set{iSubject})>0  & (yz(test_set{iSubject})>0)  );
          f2   = nnz(~y(train_set{iSubject})>0 & (yz(train_set{iSubject})>0) );
          err1 = nnz( y(test_set{iSubject})>0  ~= (yz(test_set{iSubject})>0 ));
          err2 = nnz( y(train_set{iSubject})>0 ~= (yz(train_set{iSubject})>0));
          nt1  = nnz(y(test_set{iSubject})>0);
          nt2  = nnz(y(train_set{iSubject})>0);
          nd1  = nnz(y(test_set{iSubject})<=0);
          nd2  = nnz(y(train_set{iSubject})<=0);

          if ~SMALL
            results(iii).Wz = wz;
            results(iii).Wix = uint32(ix);
            results(iii).Yz = yz;
          end

          results(iii).Wnz = uint32(wnz);
          results(iii).nvox = uint32(nv);
          results(iii).subject = uint8(iSubject);
          results(iii).cvholdout = uint8(icv);
          results(iii).finalholdout = uint8(0); % handled in parent function
          results(iii).alpha = alpha;
          results(iii).lambda = lambda;
          results(iii).diameter = 0;
          results(iii).overlap = 0;
          results(iii).shape = '';
          results(iii).nt1  = uint16(nt1);
          results(iii).nt2  = uint16(nt2);
          results(iii).nd1  = uint16(nd1);
          results(iii).nd2  = uint16(nd2);
          results(iii).h1   = uint16(h1);
          results(iii).h2   = uint16(h2);
          results(iii).f1   = uint16(f1);
          results(iii).f2   = uint16(f2);
          results(iii).err1 = uint16(err1);
          results(iii).err2 = uint16(err2);

          fprintf('cv %3d: %3d |%6.2f |%6.2f |%10d |%10d |%10.4f |%10.4f |%10d |\n', ...
            icv,iSubject,alpha_j,lambda_k,err1,err2,(h1/nt1)-(f1/nd1),(h2/nt2)-(f2/nd2),wnz);
        end
      end % lamba loop
    end % alpha loop
	end % cv loop
  fprintf('logged %d results in memory.\n', iii);
end
