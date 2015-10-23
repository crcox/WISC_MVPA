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
  addParameter(p , 'debias'         , true     );
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
  DEBIAS    = p.Results.debias;
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

  fprintf('%7s%5s %7s %7s %11s %11s %11s %11s %11s\n', '','subj','alpha','lambda','test err','train err','test diff','train diff','n vox')
  iii = 0;
  for i = 1:ncv
    icv = cvset(i);
    X = Xorig;

    test_set = cell(size(cvind));
    train_set = cell(size(cvind));
    for ii = 1:numel(cvind)
      test_set{ii}  = cvind{ii}==icv;
      train_set{ii} = ~test_set{ii};
    end

    for ii = 1:numel(X);
      switch normalize
        case 'zscore_train'
          mm = mean(X{ii}(train_set,:),1);
          ss = std(X{ii}(train_set,:),0,1);
        otherwise
          % All other cases already handled.
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
          Wz = cell(size(X));
          for ii = 1:numel(X);
            Wz{ii} = randn(size(X{ii},2),1);
          end
          info.message = 'DEBUG';
          info.iter    = 0;
        else
          switch Gtype
          case 'lasso'
            [Wz, info] = SOS_logistic( ...
              subsetAll(X, train_set), ...
              subsetAll(Y, train_set), ...
                     0, lambda(k),noG, ...
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
            opts_debias = glmnetSet(struct('alpha',0));
            for ii = 1:size(Wz)
              z = Wz{ii} ~= 0;
              if nnz(z) > 0
                cv = cvind{ii};
                cv = cv(cv~=icv);
                cv(cv>icv) = cv(cv>icv) - 1;

                debiasObj = cvglmnet(X{ii}(train_set{ii},z),Y{ii}(train_set{ii}), ...
                           'binomial',opts_debias,'class',max(cv),cv);
                opts_debias.lambda = debiasObj.lambda_min;
                debiasObj = glmnet(X{ii}(train_set{ii},z),Y{ii}(train_set{ii}),'binomial',opts_debias);
                Wz{ii}(z) = debiasObj.beta;
              end
            end
          end
        end

        k1 = sum(cellfun(@nnz,Wz))/t;

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

        for ii = 1:numel(X);
          iii = iii + 1;

          wz = Wz{ii}
          ix = find(wz);
          nv = numel(wz);
          wnz = nnz(wz);
          wz = wz(ix);

          y = Y{ii};
          yz = Yz{ii};
          h1   = nnz( y(test_set{ii})>0  & (yz(test_set{ii})>0)  );
          h2   = nnz( y(train_set{ii})>0 & (yz(train_set{ii})>0) );
          f1   = nnz(~y(test_set{ii})>0  & (yz(test_set{ii})>0)  );
          f2   = nnz(~y(train_set{ii})>0 & (yz(train_set{ii})>0) );
          err1 = nnz( y(test_set{ii})>0  ~= (yz(test_set{ii})>0 ));
          err2 = nnz( y(train_set{ii})>0 ~= (yz(train_set{ii})>0));
          nt1  = nnz(y(test_set{ii})>0);
          nt2  = nnz(y(train_set{ii})>0);
          nd1  = nnz(y(test_set{ii})<=0);
          nd2  = nnz(y(train_set{ii})<=0);

          if ~SMALL
            results(iii).Wz = wz;
            results(iii).Wix = uint32(ix);
            results(iii).Yz = Yz{ii};
          end

          results(iii).Wnz = uint32(wnz);
          results(iii).nvox = uint32(nv);
          results(iii).subject = uint8(ii);
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

          fprintf('cv %3d: %3d |%6.2f |%6.2f |%10d |%10d |%10d |%10d |%10d |\n', ...
            i,ii,alpha_j,lambda_k,err1,err2,h1-f1,h2-f2,wnz);
        end
      end % lamba loop
    end % alpha loop
	end % cv loop
  fprintf('logged %d results in memory.\n', iii);
end
