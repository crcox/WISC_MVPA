% Run jobs in this directory
addpath('C:\Users\mbmhscc4\MATLAB\src\WholeBrain_MVPA\src');
addpath('C:\Users\mbmhscc4\MATLAB\src\WholeBrain_MVPA\dependencies\jsonlab\');
dlist = selectbyfield(dir('.'), 'isdir', 1);
z = ~cellfun('isempty', regexp({dlist.name}, '[0-9]+'));
job_dirs = dlist(z);

for i = 1:numel(job_dirs)
    jdir = job_dirs(i).name;
    cd(jdir);
    WholeBrain_MVPA;
    cd('..');
end