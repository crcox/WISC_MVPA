function test_soglasso()
  format compact
  warning('off');

  rng(10);

  tic;
  fprintf('Testing the SOGlasso code \n')

  n = 800;
  t = 5;
  k = 5;
  l = 30;

  fprintf('forming overlapping groups...')
  L = 50;B = 50;
  G = cell(0);
  for g = 1:L
      if g==1
          G{g} = (g-1)*B+1:g*B;
          lastind = g*B;
      else
          G{g} = lastind-5:lastind-5+B-1;
          lastind = lastind-2+B-1;
      end
  end
  fprintf('done. time = %f \n',toc);
  p = lastind;


  fprintf('signal dimension                  = %d \n',p);
  fprintf('number of measurements            = %d \n',n);
  fprintf('number of tasks                   = %d \n',t);
  fprintf('number of groups                  = %d \n',L);
  fprintf('number of active groups           = %d \n',k);
  fprintf('number of active entries per group= %d \n',l);



  A = randn(n,p);
  [A_tilde, groups, group_arr] = replicate_matrix(A,G);
  fprintf('Features replicated. Time = %f \n',toc);

  W = sparse(zeros(p,t));
  acts = randsample(L,k);
  for g = 1:length(acts)
      grp = G{acts(g)};
      temp = zeros(B*t,1);
      actinds = randsample(B*t,l);
      temp(actinds) = 2*rand(l,1)-1;
      temp = reshape(temp,B,t);
      W(grp,:) = W(grp,:) + temp;
  end

  fprintf('\n _______ LEAST SQUARES LOSS ________ \n')
  y = A*W + 0.02*randn(n,t);
  for task = 1:t
      Y{task} = y(:,task);
      X{task} = A_tilde;
  end



  alpha = 0.1; lambda = 0.25;
  [Wtemp,obj] = SOG_least_squares(X, Y,alpha,lambda,groups,group_arr);

  % debias
  Wopt = zeros(size(Wtemp));
  for task = 1:t
      inds = find(Wtemp(:,task)~=0);
      x = X{task};
      x = x(:,inds);
      Wopt(inds,task) = x\Y{task};
  end

  What = sparse(zeros(size(W)));
  for g = 1:L
      grpnew = find(groups==g);
      grpold = G{g};
      What(grpold,:) = Wopt(grpnew,:) + What(grpold,:);
  end

  error_ls = norm(What-W,'fro')^2/numel(W);
  fprintf('LS error = %f . Time = %f \n',error_ls,toc);


  fprintf('\n _______________ LOGISTIC LOSS ________________ \n')
  for task = 1:t
      Y{task} = sign(Y{task});
  end

  alpha = 0.1; lambda = 0.25;
  opts.maxiter = 100;
  [Wtemp,obj] = SOG_logistic(X, Y,alpha,lambda,groups,group_arr,'opts',opts);

  % debias
  Wopt = zeros(size(Wtemp));
  for task = 1:t
      inds = find(Wtemp(:,task)~=0);
      x = X{task};
      x = x(:,inds);
      [ Wopt(inds,task),obj] =...
          SOG_logistic({x}, Y(task),0,lambda,groups,group_arr);
  end

  What = sparse(zeros(size(W)));
  for g = 1:L
      grpnew = find(groups==g);
      grpold = G{g};
      What(grpold,:) = Wopt(grpnew,:) + What(grpold,:);
  end

  error_log = norm(What-W,'fro')^2/numel(W);
  fprintf('LOGISTIC error = %f . Time = %f \n',error_log,toc);
end


