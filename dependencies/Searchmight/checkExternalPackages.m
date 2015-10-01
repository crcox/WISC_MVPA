function [] = checkExternalPackages()

% test for all our builds (doesn't test specifically for the architecture in *this* machine,
% still have to figure out how MATLAB would do that)
[location_mexglx]  = which('svmtrain.mexglx');
[location_mexmaci] = which('svmtrain.mexmaci');
[location_mexa64]  = which('svmtrain.mexa64');

test_svmtrain = isempty(location_mexglx) & isempty(location_mexmaci) & isempty(location_mexa64);

if test_svmtrain; fprintf('ERROR: cannot find a build of libsvm svmtrain\n');return; end

[location_mexglx]  = which('svmpredict.mexglx');
[location_mexmaci] = which('svmpredict.mexmaci');
[location_mexa64]  = which('svmpredict.mexa64');

test_svmpredict = isempty(location_mexglx) & isempty(location_mexmaci) & isempty(location_mexa64);

if test_svmpredict; fprintf('ERROR: cannot find a build of libsvm svmpredict\n');return; end

fprintf('OK: all external packages are there!\n');
