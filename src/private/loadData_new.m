function X = loadData_new(datafiles, data_var, FMT_subjid)
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
    X = struct('filename',datafiles,'subject',[],data_var,[]);
    N = length(datafiles);
    for i = 1:N
        x = datafiles{i};
        fprintf('Loading %s from  %s...\n', data_var, x);
        tmp = load(datafiles{i}, data_var);
        X(i).(data_var) = tmp.(data_var); clear tmp;
        X(i).subject = extractSubjectID(x, FMT_subjid);
        if hasWindowInfo(x)
            s = extractWindowInfo(x);
            X(i).BoxCar = s.BoxCar;
            X(i).WindowStart = s.WindowStart;
            X(i).WindowSize = s.WindowSize;
        end
    end
end

function s = extractWindowInfo(x)
    a = strfind(x,'BoxCar');
    y = sscanf(x(a:end), 'BoxCar/%d/WindowStart/%d/WindowSize/%d');
    s = struct('BoxCar',y(1),'WindowStart',y(2),'WindowSize',y(3));
end

function b = hasWindowInfo(x)
    y = {'BoxCar','WindowSize','WindowStart'};
    b = all(ismember(y,strsplit(x, '/'))) || all(ismember(y,strsplit(x, '\')));
end