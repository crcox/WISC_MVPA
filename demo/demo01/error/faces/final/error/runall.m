addpath('~/src/WholeBrain_MVPA/src')
addpath('~/src/WholeBrain_MVPA/dependencies/jsonlab/')
addpath('~/src/WholeBrain_RSA/util')
for i = 1:9
  cd(num2str(i));
  WholeBrain_MVPA
  cd('..');
end
  