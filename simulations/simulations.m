%%
SOSLASSODIR = '~/src/soslasso';
if ~exist(SOSLASSODIR,'dir')
  SOSLASSODIR = 'C:\Users\chris\Documents\soslasso';
end
addpath(fullfile(SOSLASSODIR,'src'));
addpath(fullfile(SOSLASSODIR,'util'));
addpath(fullfile(SOSLASSODIR,'simulations'));

%%
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
y= cell2mat(Y');
%%
group_arr = ndmovingwindow(nvox, 'size', 10);
groups = repmat((1:size(group_arr,1))',1,size(group_arr,2));
groups = groups(:)';
G = repmat(mat2cell(group_arr,ones(size(group_arr,1),1),size(group_arr,2)),1,numel(X));
GroupInfo.groups = groups;
GroupInfo.group_arr = group_arr;
GroupInfo.RepIndex = 1:nvox;
%%
CV = repmat({bsxfun(@eq, 1+mod(1:nitem,4)', 1:4)},nsubj,1);
lambda = linspace(2,5,10);
alpha = 0.01;
%opts = struct(); % defaults will be set internally.
cv = cellfun(@(x) x(:,1), CV, 'Unif', 0);
Bz = cell(1,length(lambda));
ErrorC = cell(1,numel(lambda));
nzv = zeros(numel(lambda),numel(X));
for i = 1:length(lambda)
  disp(i)
  Bz{i} = SOS_logistic( ...
  subsetAll(X,[],[]), ...
  subsetAll(Y,[],[]), ...
  alpha, ...
  lambda(i), ...
  G);
end
Yz = cell(1,numel(Bz));
for j = 1:numel(Bz)
  disp(j)
  Yz{j} = zeros(nitem,nsubj);
  bz = Bz{j};
  for i = 1:nsubj
    z = bz{i}~=0;
    if nnz(z)
      cvopts = glmnetSet();
      cvopts.alpha = 0;
      cvopts.intr = 0;
      cvobj = cvglmnet(X{i}(:,z),Y{i},'binomial',cvopts,'class',10,repmat(1:10,1,10)');
      opts = glmnetSet();
      opts.lambda = cvobj.lambda_min;
      opts.alpha = 0;
      opts.intr = 0;
      fitobj = glmnet(X{i}(:,z),Y{i},'binomial',opts);
      bz{i}(z) = fitobj.beta;
    end
    Yz{j}(:,i) = X{i} * bz{i};
    nzv(j,i) = nnz(bz{i});
  end
  ErrorC{j} = mean(y~=(Yz{j}>0))';
end
subplot(2,1,1)
plot(cell2mat(ErrorC)')
subplot(2,1,2)
plot(nzv)