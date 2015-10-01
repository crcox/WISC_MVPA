%
% set up paths for this toolbox and its dependences
%

function [] = setupPathsSearchmightToolbox()

current = pwd;
tmp = which('setupPathsSearchmightToolbox');
SearchmightToolboxLocation = fileparts(tmp);
cd(SearchmightToolboxLocation);

% get bearings
[status,locationTxt]  = unix('pwd');
%[locationTxt,discard] = strtok(locationTxt,' ');
locationTxt = sprintf('%s',locationTxt(1:(end-1)));
[status,systemTxt]    = unix('uname -s');
[status,processorTxt] = unix('uname -p');
[status,releaseTxt]   = unix('uname -r');

% add architecture specific packages
fprintf('+%s+\n',locationTxt);
targetDirectory = sprintf('ExternalPackages.%s_%s',systemTxt(1:(end-1)),processorTxt(1:(end-1)));
dtxt = sprintf('%s/CoreToolbox/%s/libsvm',locationTxt,targetDirectory);

fprintf('+%s+\n',dtxt);
addpath(dtxt);

cd(current);
