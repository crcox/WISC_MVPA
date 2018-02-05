function X = loadData_new(datafiles, data_var)
% LOADDATA Load example-by-feature matrices for multiple subjects.
%
% A simple function that allows the matrices to easily be loaded into a
% arbitrarily named variable as a cell array in the workspace.
%
% datafiles : A cell array of paths to .mat files to load.
% data_var  : The name of the specific variable from each .mat file. Assumes
%             that variable will have the same name across datafiles.
%
% Chris Cox 24/08/2017
    datafiles  = ascell(datafiles);
    N         = length(datafiles);
    X         = cell(N,1);
    for i = 1:N
        fprintf('Loading %s from  %s...\n', data_var, datafiles{i});
        tmp       = load(datafiles{i}, data_var);
        X{i}      = tmp.(data_var); clear tmp;
    end
end
