function X = loadData_new(datafiles, data_var, metafile, metadata_varname, FMT_subjid, filter_labels, target_label, target_type, sim_source, sim_metric, cvscheme, finalholdoutInd)
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
    X = struct('filename',datafiles,'subject',[],data_var,[],'rowfilters',[],'colfilters',[]);
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
        M = selectbyfield(metadata, 'subject', X(i).subject);
        if numel(M) == 1
            T = selectbyfield(M.targets, ...
                'label', target_label, ...
                'type', target_type, ...
                'sim_source', sim_source, ...
                'sim_metric', sim_metric);
            X(i).targets = T;
            
            RF = selectbyfield(M.filters, 'label', filter_labels, 'dimension', 1);
            if strcmpi('notes', fieldnames(RF))
                FHO = struct( ...
                    'label', 'finalholdout', ...
                    'orientation', 1, ...
                    'filter', M.cvind(:,cvscheme) == finalholdoutInd, ...
                    'note', 'autogen');
            else
                FHO = struct( ...
                    'label', 'finalholdout', ...
                    'dimension', 1, ...
                    'filter', M.cvind(:,cvscheme) == finalholdoutInd);
            end
            if any(strcmpi('finalholdout',{RF.label}))
                RF = replacebyfield(RF,FHO,'label','finalholdout','dimension',1);
            else
                RF(end+1) = FHO; %#ok<AGROW> (it's not actually growing in a loop)
            end
            CF = selectbyfield(M.filters, 'label', filter_labels, 'dimension', 2);
            for j = 1:numel(RF), RF(j).filter = RF(j).filter(:)'; end
            for j = 1:numel(CF), CF(j).filter = CF(j).filter(:)'; end
            X(i).rowfilters = RF;
            X(i).colfilters = CF;
            
            X(i).coords = M.coords;

        else
            error('WholeBrain_MVPA:loadData:SubjectSelect','Subject ID %s does not uniquely match a metadata entry.', num2cell(X(i).subject));
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