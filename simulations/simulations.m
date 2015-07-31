addpath('~/src/soslasso/src');
nvox = 1000;
nsubj = 10;
nitem = 100;

b = randn(nvox,nsubj);
b(rand(nvox,nsubj)<0.95) = 0;
b = mat2cell(b, nvox, ones(1,nsubj));
cellfun(@nnz, b);

X = mat2cell(randn(nsubj*nitem,nvox),repmat(nitem,nsubj,1),nvox);

Y = cell(nsubj,1);
for i = 1:nsubj
  Y{i} = (X{i}*b{i})>0;
end

group_arr = ndmovingwindow(nvox, 'size', 10);
groups = repmat((1:size(group_arr,1))',1,size(group_arr,2));
groups = groups(:)';
GroupInfo.G = mat2cell(group_arr,ones(size(group_arr,1),1),size(group_arr,2));
GroupInfo.groups = groups;
GroupInfo.group_arr = group_arr;
GroupInfo.RepIndex = 1:nvox;

CV = repmat({bsxfun(@eq, 1+mod(1:nitem,4)', 1:4)},nsubj,1);

lambda = linspace(0,0.03,100);
alpha = 0.1;
opts = struct(); % defaults will be set internally.
cv = cellfun(@(x) x(:,1), CV, 'Unif', 0);
Bz = cell(1,length(lambda));
err = @(yz) mean(mean(y~=(yz>0)));
errc = @(yz) mean(y~=(yz>0));
for i = 1:length(lambda)
  Bz{i} = SOG_logistic( ...
  subsetAll(X,[],[]), ...
  subsetAll(Y,[],[]), ...
  alpha, ...
  lambda(i), ...
  groups, ...
  group_arr, ...
  opts);
end
Yz = cell(1,numel(Bz));
for j = 1:numel(Bz)
  Yz{j} = zeros(nitem,nsubj);
  bz = Bz{j};
  for i = 1:nsubj
    Yz{j}(:,i) = X{i} * bz(:,i);
  end
end
ErrorC = cellfun(errc, Yz, 'Unif', 0);
Error = cellfun(err, Yz);
hold on
plot(Error)
%% cvsoslasso
lambda = linspace(0,100,100);
for i = 1:length(lambda)
  Bz2{i} = overlap_2stage(Y,X,GroupInfo,lambda(i),1,init_opts(struct()));
end

Yz = cell(1,numel(Bz2));
for j = 1:numel(Bz2)
  Yz{j} = zeros(nitem,nsubj);
  bz = Bz2{j};
  for i = 1:nsubj
    Yz{j}(:,i) = X{i} * bz(:,i);
  end
end
ErrorC = cellfun(errc, Yz, 'Unif', 0);
plot(cell2mat(ErrorC'))
y = cell2mat(Y');
yz = cell2mat(Yz);
mean(mean(y~=(yz>0)))
sum(bz~=0)

b = cell2mat(b);
subplot(1,2,1)
imagesc(b)
subplot(1,2,2)
imagesc(bz)