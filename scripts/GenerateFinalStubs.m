% Setup an analysis based on the parameter search

% This template should be copied into the directory structure of a specific
% analysis, and modified to fit needs.

addpath('C:\Users\mbmhscc4\MATLAB\src\WholeBrain_MVPA\util');
addpath('C:\Users\mbmhscc4\MATLAB\src\WholeBrain_MVPA\dependencies\jsonlab\');
addpath(genpath('C:\Users\mbmhscc4\MATLAB\Toolboxes\YAMLMatlab_0.4.3'))
% Note: I am using a modified version of WriteYaml, where scan numeric is
% defined as:
% function result = scan_numeric(r)
%     if isempty(r)
%         result = java.util.ArrayList();
%     elseif isinteger(r) % CRC ADDITION
%         result = java.lang.Integer(r);
%     else
%         result = java.lang.Double(r);
%     end
% end

%% Load the results of the parameter search
[~,tune,n] = HTCondorLoad('HB_0');

%% Identify the best configurations
MODALITY = 'visual';

T = struct2table(tune);
nrow = size(T,1);
T.modality = repmat({MODALITY},nrow,1);
HTCondor_struct2csv(tune, 'HB_0.csv', 'Skip', {'Uz','Uix','nz_rows','Sz','Cz','coords','jobdir'});

BC = BestCfgByCondition(T,{'lambda','lambda1'},'err1',{'nzvox'},{'subject','modality'});
writetable(BC,'BestCfgByCondition.csv')

%% Read the config file for the parameter search ("tune stub")
tune_stub = ReadYaml('HB_0\stub.yaml');
[fpath,~,~] = fileparts(tune_stub.data{1});
fformat = 's%02d_sessions.mat'; % tune_stub.subject_id_fmt;
subject_data_files = arrayfun(@(i) sprintf(fformat, i), BC.subject', 'unif', false);

%% Based on the tune stub and best config, build the final stub
final_stub = tune_stub;
final_stub.bias = int32(final_stub.bias);
final_stub.lambda = BC.lambda';
final_stub.lambda1 = BC.lambda1';
final_stub.data = fullfile(fpath, subject_data_files);
final_stub.data = strrep(final_stub.data, '\','/');
final_stub.cvscheme = int32(1);
final_stub.cvholdout = int32(0);
final_stub.finalholdout = int32(0);
final_stub.SmallFootprint = false;
final_stub.EXPAND = {{'lambda', 'lambda1', 'data'}};
final_stub = rmfield(final_stub, 'HYPERBAND');

WriteYaml('stub_final.yaml', final_stub);

%% Based on the final stub, build the permutations stub
permfile_name = 'PERMUTATION_STRUCTURE_FULL.mat';
perm_stub = final_stub;
perm_stub.PermutationTest = true;
perm_stub.PermutationMethod = 'manual';
perm_stub.PermutationIndex = fullfile(fpath, permfile_name);
perm_stub.PermutationIndex = strrep(perm_stub.PermutationIndex, '\','/');
perm_stub.perm_varname = 'PERMUTATION_INDEX';
perm_stub.RandomSeed = reshape(int32(1:100),10,10)';
perm_stub.EXPAND = {{'lambda', 'lambda1', 'data'}, 'RandomSeed'};
perm_stub.URLS = [perm_stub.URLS,{'PermutationIndex'}];

WriteYaml('stub_permutations.yaml', perm_stub);

%% Generate stub files that can be used locally
fpath_local = 'D:\MRI\SoundPicture\data\MAT\avg\bysession';
[~, meta_name, ~] = fileparts(final_stub.metadata);
final_stub_local = final_stub;
final_stub_local.data = fullfile(fpath_local, subject_data_files);
final_stub_local.metadata = fullfile(fpath_local, meta_name);
final_stub_local = rmfield(final_stub_local, {'executable','wrapper','COPY','URLS'});
WriteYaml('stub_final_local.yaml', final_stub_local);

perm_stub_local = perm_stub;
perm_stub_local.data = fullfile(fpath_local, subject_data_files);
perm_stub_local.metadata = fullfile(fpath_local, meta_name);
perm_stub_local.PermutationIndex = fullfile(fpath_local, permfile_name);
perm_stub_local = rmfield(perm_stub_local, {'executable','wrapper','COPY','URLS'});
WriteYaml('stub_permutations_local.yaml', perm_stub_local);

