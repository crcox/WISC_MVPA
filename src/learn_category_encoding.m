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
  addParameter(p , 'verbose'        , false    );
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
  VERBOSE   = p.Results.verbose;
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
  results.confusion1 = [];
  results.confusion2 = [];
  results.iterations = struct();

  % Preallocate
  results(t*ncv*nlam*nalpha).Wz = [];

  % Permute if requested
  fprintf('PermutationTest: %d\n', PermutationTest);
  if PermutationTest
    if iscell(Y)
      for i = 1:numel(Y)
        for ic = unique(cvind{i})'
          fprintf('Permuting CV %d...\n', ic);
          y = Y{i}(cvind{i}==ic);
          n = size(y,1);
          if VERBOSE
            fprintf('Permuting %d rows of C, independently by its %d columns.\n', n, r);
            fprintf('First 10 rows of C, before shuffling.\n')
            disp(y)
          end
          permix = randperm(n);
          Y{i}(cvind{i}==ic) = y(permix, :);
          if VERBOSE
            fprintf('First 10 rows of C, after shuffling.\n')
            disp(y(permix))
          end
        end
      end
    else
      for ic = unique(cvind)'
        fprintf('Permuting CV %d...\n', ic);
        y = Y(cvind==ic);
        n = size(y,1);
        if VERBOSE
          fprintf('Permuting %d rows of Y.\n', n, r);
          fprintf('First 10 rows of Y, before shuffling.\n')
          disp(y)
        end
        permix = randperm(n);
        Y(cvind==ic) = y(permix, :);
        if VERBOSE
          fprintf('First 10 rows of Y, after shuffling.\n')
          disp(y(permix))
        end
      end
    end
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
          mm1 = mean(X{iSubject}(test_set{iSubject},:),1);
          ss1 = std(X{iSubject}(test_set{iSubject},:),0,1);
          mm2 = mean(X{iSubject}(train_set{iSubject},:),1);
          ss2 = std(X{iSubject}(train_set{iSubject},:),0,1);
        otherwise
          % All other cases already handled.
          [mm1,mm2] = deal(0);
          [ss1,ss2] = deal(1);
      end
      X{iSubject}(test_set{iSubject},:) = bsxfun(@minus,X{iSubject}(test_set{iSubject},:), mm1);
      X{iSubject}(test_set{iSubject},:) = bsxfun(@rdivide,X{iSubject}(test_set{iSubject},:), ss1);
      X{iSubject}(train_set{iSubject},:) = bsxfun(@minus,X{iSubject}(train_set{iSubject},:), mm2);
      X{iSubject}(train_set{iSubject},:) = bsxfun(@rdivide,X{iSubject}(train_set{iSubject},:), ss2);
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
          case 'smlr'
            [Wz, info] = smlr_mvpatoolbox( ...
              subsetAll(X, train_set), ...
              subsetAll(Y, train_set), ...
              lambda,...
              'maxiter' , 1e4, ...
              'tol'     , 1e-5, ...
              'verbose' , true, ...
              'W0'      ,   []);

          case 'lasso_glmnet'
            if iscell(Y)
              for iy = 1:numel(Y)
                Y{iy} = squeeze_index(Y{iy});
              end
            else
              Y = squeeze_index(Y);
            end
            [Wz, info] = lasso_glmnet( ...
              subsetAll(X, train_set), ...
              subsetAll(Y, train_set), ...
                       1, lambda,...
              'cvind'   , subsetAll(cvind, train_set),...
              'maxiter' , 1000, ...
              'tol'     , 1e-8, ...
              'W0'      ,   []);

          case 'iterlasso_glmnet'
            % Ensure that y indexes are consecutive, with min=1
            if iscell(Y)
              for iy = 1:numel(Y)
                Y{iy} = squeeze_index(Y{iy});
              end
            else
              Y = squeeze_index(Y);
            end
            [Wz, info, Iterations] = iterlasso_glmnet( ...
              subsetAll(X, train_set), ...
              subsetAll(X, test_set), ...
              subsetAll(Y, train_set), ...
              subsetAll(Y, test_set), ...
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
          nz_rows = any(wz,2);
          ix = find(nz_rows);
          nv = size(wz,1);
          wnz = nnz(nz_rows);
          wz = wz(nz_rows,:);
          
          y = Y{iSubject};
          ytest = y(test_set{iSubject});
          ytrain = y(train_set{iSubject});
          
          cinds = unique(y);
          cinds = cinds(:)'; % force to row vec
          m = numel(cinds);
          if m > 2;
            isMultinomial = 1;
            ytrain_m = bsxfun(@eq,ytrain(:),cinds);
            ytest_m = bsxfun(@eq,ytest(:),cinds);
          else
            isMultinomial = 0;
            ytrain_m = ytrain;
            ytest_m = ytest;
          end

          yz_m = X{iSubject} * Wz{iSubject};
          yztest_m = yz_m(test_set{iSubject},:);
          yztrain_m = yz_m(train_set{iSubject},:);
          
          h1   = sum((ytest_m > 0) & (yztest_m > 0));
          f1   = sum(~(ytest_m > 0) & (yztest_m > 0));
          err1 = sum((ytest_m > 0) ~= (yztest_m > 0));
          nt1  = sum(ytest_m > 0);
          nd1  = sum(ytest_m <= 0);

          h2   = sum((ytrain_m > 0) & (yztrain_m > 0));
          f2   = sum(~(ytrain_m > 0) & (yztrain_m > 0));
          err2 = sum((ytrain_m > 0) ~= (yztrain_m > 0));
          nt2  = sum(ytrain_m > 0);
          nd2  = sum(ytrain_m <= 0);

          if isMultinomial
            [~,yztest] = max(yztest_m,[],2);
            [~,yztrain] = max(yztrain_m,[],2);
            confusion1 = confusionmat(double(ytest),yztest);
            confusion2 = confusionmat(double(ytrain),yztrain);
          else
            [confusion1,confusion2] = deal([]);
          end
          
          
%           if isMultinomial
%               [~,yz] = max(yz,[],2);
%               err1 = nnz(y(test_set{iSubject})  ~= yz(test_set{iSubject}));
%               err2 = nnz(y(train_set{iSubject}) ~= yz(train_set{iSubject}));
%               confusion1 = confusionmat(y(test_set{iSubject}),yz(test_set{iSubject}));
%               confusion2 = confusionmat(y(train_set{iSubject}),yz(train_set{iSubject}));
%               [h1,h2,f1,f2,nt1,nt2,nd1,nd2] = deal([]);
%           else
%               h1   = nnz( y(test_set{iSubject})>0  & (yz(test_set{iSubject})>0)  );
%               h2   = nnz( y(train_set{iSubject})>0 & (yz(train_set{iSubject})>0) );
%               f1   = nnz(~y(test_set{iSubject})>0  & (yz(test_set{iSubject})>0)  );
%               f2   = nnz(~y(train_set{iSubject})>0 & (yz(train_set{iSubject})>0) );
%               err1 = nnz( y(test_set{iSubject})>0  ~= (yz(test_set{iSubject})>0 ));
%               err2 = nnz( y(train_set{iSubject})>0 ~= (yz(train_set{iSubject})>0));
%               nt1  = nnz(y(test_set{iSubject})>0);
%               nt2  = nnz(y(train_set{iSubject})>0);
%               nd1  = nnz(y(test_set{iSubject})<=0);
%               nd2  = nnz(y(train_set{iSubject})<=0);
%               [confusion1,confusion2] = deal([]);
%           end

          if ~SMALL
            results(iii).Wz = wz;
            results(iii).Wix = uint32(ix);
            results(iii).Yz = yz_m;
          end

          results(iii).Wnz = uint32(wnz);
          results(iii).nvox = uint32(nv);
          results(iii).subject = uint8(iSubject);
          results(iii).cvholdout = uint8(icv);
          results(iii).finalholdout = uint8(0); % handled in parent function
          results(iii).alpha = alpha;
          results(iii).lambda = lambda;
          results(iii).diameter = 0; % handled in parent function
          results(iii).overlap = 0; % handled in parent function
          results(iii).shape = ''; % handled in parent function
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
          results(iii).confusion1 = confusion1;
          results(iii).confusion2 = confusion2;
          switch lower(regularization)
            case 'iterlasso_glmnet'
              results(iii).iterations = Iterations;
          end

          if isMultinomial
            fprintf('cv %3d: %3d |%6.2f |%6.2f |%10.4f |%10.4f |%10.4f |%10.4f |%10d |\n', ...
              icv,iSubject,alpha_j,lambda_k,mean(err1),mean(err2),mean((h1/nt1)-(f1/nd1)),mean((h2/nt2)-(f2/nd2)),wnz);
          else
            fprintf('cv %3d: %3d |%6.2f |%6.2f |%10d |%10d |%10.4f |%10.4f |%10d |\n', ...
              icv,iSubject,alpha_j,lambda_k,err1,err2,(h1/nt1)-(f1/nd1),(h2/nt2)-(f2/nd2),wnz);
          end
        end
      end % lamba loop
    end % alpha loop
	end % cv loop
  fprintf('logged %d results in memory.\n', iii);
end
function [s] = squeeze_index(x)
  u = unique(sort(x));
  s = zeros(size(x));
  for i = 1:numel(u);
      s(x==u(i)) = i;
  end
end