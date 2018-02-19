function subjectArray = loadData_new(datafiles, data_var, metafile, metadata_varname, FMT_subjid, filter_labels, target_label, target_type, sim_source, sim_metric, cvscheme, finalholdoutInd)
% LOADDATA Load example-by-feature matrices for multiple subjects. 
%
% A simple function that allows the matrices to easily be loaded into a
% arbitrarily named variable as a cell array in the workspace.
%
% datafiles : A cell array of paths to .mat files to load.
% data_var  : The name of the specific variable from each .mat file. Assumes
%             that variable will have the same name across datafiles.
%
% Chris Cox 18/02/2018
    StagingContainer = load(metafile, metadata_varname);
    metadata = StagingContainer.(metadata_varname); clear StagingContainer;
%     subject_label = {metadata.subject};
    [metadata, ~] = subsetMetadata(metadata, datafiles, FMT_subjid);

%     tmpS = selectTargets(metadata, target_type, target_label, sim_source, sim_metric, rowfilter(subjix));
%     S = cell(max(subjix),1);
%     for i = 1:numel(tmpS);
%         S(subjix(i)) = tmpS(i);
%     end
%     clear tmpS;
    
    datafiles  = ascell(datafiles);
    subjectArray = repmat(Subject, 1, numel(datafiles));
    
    N = length(datafiles);
    
    for i = 1:N
        x = datafiles{i};
        fprintf('Loading %s from  %s...\n', data_var, x);
        X = struct( ...
            'filename',datafiles{i}, ... % set right away
            'subject',[], ...
            'label',data_var, ... % set right away
            'data',[], ...
            'cvscheme', [], ...
            'rowfilters',[], ...
            'colfilters',[], ...
            'targets', [], ...
            'coords',[]);
        tmp = load(X.filename, X.label);
        X.data = tmp.(data_var); clear tmp;
        X.subject = extractSubjectID(x, FMT_subjid);
        if hasWindowInfo(x)
            s = extractWindowInfo(x);
            X.BoxCar = s.BoxCar;
            X.WindowStart = s.WindowStart;
            X.WindowSize = s.WindowSize;
        end
        M = selectbyfield(metadata, 'subject', X.subject);
        X.cvscheme = M.cvind(:,cvscheme);
        if numel(M) == 1
            T = selectbyfield(M.targets, ...
                'label', target_label, ...
                'type', target_type, ...
                'sim_source', sim_source, ...
                'sim_metric', sim_metric);
            X.targets = T;
            
            RF = selectbyfield(M.filters, 'label', filter_labels, 'dimension', 1);
            if strcmpi('notes', fieldnames(RF))
                FHO = struct( ...
                    'label', 'finalholdout', ...
                    'orientation', 1, ...
                    'filter', X.cvscheme ~= finalholdoutInd, ...
                    'note', 'autogen');
            else
                FHO = struct( ...
                    'label', 'finalholdout', ...
                    'dimension', 1, ...
                    'filter', X.cvscheme ~= finalholdoutInd);
            end
            if any(strcmpi('finalholdout',{RF.label}))
                RF = replacebyfield(RF,FHO,'label','finalholdout','dimension',1);
            else
                RF(end+1) = FHO; %#ok<AGROW> (it's not actually growing in a loop)
            end
            CF = selectbyfield(M.filters, 'label', filter_labels, 'dimension', 2);
            for j = 1:numel(RF), RF(j).filter = RF(j).filter(:)'; end
            for j = 1:numel(CF), CF(j).filter = CF(j).filter(:)'; end
            X.rowfilters = RF;
            X.colfilters = CF;
            X.coords = M.coords;
            subjectArray(i) = Subject(X);

        else
            error('WholeBrain_MVPA:loadData:SubjectSelect','Subject ID %s does not uniquely match a metadata entry.', num2cell(X.subject));
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