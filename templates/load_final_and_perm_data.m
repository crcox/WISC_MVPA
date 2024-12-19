function load_final_and_perm_data(path)
    % the path here should be the project directory - it should CONTAIN
    % (not BE) a directory called \final or \perm.
    addpath('C:\Users\slfri\ownCloud\7T_WISC_MVPA\scripts\WISC_MVPA\src')
    addpath('C:\Users\slfri\ownCloud\7T_WISC_MVPA\scripts\WISC_MVPA\util')
    addpath('C:\Users\slfri\ownCloud\7T_WISC_MVPA\scripts\WISC_MVPA\dependencies\jsonlab')
    addpath('C:\Users\slfri\ownCloud\7T_WISC_MVPA\scripts') 

    % skip any matrices that you do not need to write to .csv (these are still
    % written to .mat). These include coefficients (Uz), predicted coordinates
    % (Cz), indices of nonzero voxels (nz_rows), coordinates of nonzero voxels
    % (coords) and perhaps filters (filters). 
    % N.B. while it might seem nice to have coefficients and coordinates saved
    % to .csv, the way that they are nested in the MATLAB table makes this very
    % difficult and the output illegible and odd - so inspect these within
    % MATLAB if needed. 
    if ~isempty(strfind(path,'correlation'))
        SKIP = {'Uz', 'Cz', 'nz_rows', 'coords'};
    elseif ~isempty(strfind(path,'classification'))
        SKIP = {'Wz', 'Yz', 'nz_rows', 'coords'};
    end

    % load the final directory tree
    Tallcv = load_from_condor([path,'\final\']);
    % outputs of SOSLASSO have some of the character variables constructed
    % in a weird way. Correct this
    if ~isempty(strfind(path,'classification'))
        for i = 1:size(Tallcv,1)
            tmp = table2array(Tallcv(i,7));
            Tallcv(i,5) = {tmp{1,1}'};
            tmp = table2array(Tallcv(i,8));
            Tallcv(i,6) = {tmp{1,1}'};
        end
    end
    % save the output as a .mat file (including coefficients, Uz, and predicted
    % coordinates, Cz)
    save([path,'\final_performance.mat'],'Tallcv','-v7.3');
    % save the abridged output as a .csv file
    writetable(removevars(Tallcv, SKIP),[path,'\final_performance.csv']);
    
    % load the perm directory tree
    Tallcv = load_from_condor([path,'\perm\']);
    % outputs of SOSLASSO have some of the character variables constructed
    % in a weird way. Correct this
    if ~isempty(strfind(path,'classification'))
        for i = 1:size(Tallcv,1)
            tmp = table2array(Tallcv(i,7));
            Tallcv(i,5) = {tmp{1,1}'};
            tmp = table2array(Tallcv(i,8));
            Tallcv(i,6) = {tmp{1,1}'};
        end
    end
    % save the output as a .mat file (including coefficients, Uz, and predicted
    % coordinates, Cz)
    save([path,'\perm_performance.mat'],'Tallcv','-v7.3');
    % save the abridged output as a .csv file
    writetable(removevars(Tallcv, SKIP),[path,'\perm_performance.csv']);
end