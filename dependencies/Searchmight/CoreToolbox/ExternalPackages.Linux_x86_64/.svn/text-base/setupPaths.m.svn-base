[status,locationTxt]  = unix('pwd'); locationTxt = locationTxt(1:(end-1));
[status,systemTxt]    = unix('uname -s');
[status,processorTxt] = unix('uname -p');
[status,releaseTxt]   = unix('uname -r');

% set up our own path
addpath(locationTxt);

%sprintf('addpath %s/cvx/builtins',locationTxt)
%return;
eval(sprintf('addpath %s/cvx/builtins',locationTxt));
eval(sprintf('addpath %s/cvx/commands',locationTxt));
eval(sprintf('addpath %s/cvx/functions',locationTxt));
eval(sprintf('addpath %s/cvx/lib',locationTxt));
eval(sprintf('addpath %s/cvx/structures',locationTxt));
eval(sprintf('addpath %s/libsvm',locationTxt));

