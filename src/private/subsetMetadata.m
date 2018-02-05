function [metadata,subjix] = subsetMetadata(metadata,datafiles,FMT_subjid)
% SUBSETMETADATA Extract subject ID from filename and match METADATA.subject
%
% Metadata for all subjects in a project is stored in a single structured
% array. Each sub-array contains a field called 'subject', which contains
% either an integer or string that uniquely identifies the subject and
% systematically matches a portion of the filename that references the
% example-by-voxel data matrix for that subject.
%
% Example 1:
% ---------
% metadata = struct('subject', {100,101}, 'targets', struct());
% datafile = {'path/to/DATA101_avg.mat'};
% fmt = 'DATA%03d_avg.mat';
%
% [M, subjix] = subsetMetadata(metadata, datafile, fmt);
% M =
%   subject : 101
%
% subjix =
%   2
%
% Example 2:
% ---------
% metadata = struct('subject', {'DATA100','DATA101'}, 'targets', struct());
% datafile = {'path/to/DATA100_avg.mat'};
% % Because subject ID is a string, FMT_subjid is not needed.
%
% [M, subjix] = subsetMetadata(metadata, datafile);
% M =
%   subject : 'DATA100'
%
% subjix =
%   1
%
% Chris Cox 24/08/2017
    datafiles = ascell(datafiles);
    N         = length(datafiles);
    subjix    = zeros(1,N);
    for i = 1:N
        if ischar(metadata(1).subject)
            [~,subjcode,~] = fileparts(datafiles{i});
            z = strcmp({metadata.subject}, subjcode);
            if ~any(z)
                z = false(1,numel(metadata));
                for j = 1:numel(metadata)
                    n = length(metadata(j).subject);
                    z(j) = strncmp(metadata(j).subject, subjcode, n);
                end
            end
            subjix(i) = find(z);
        else
            subjid    = extractSubjectID(datafiles{i}, FMT_subjid);
            subjix(i) = find([metadata.subject] == subjid);
        end
    end
    metadata = metadata(subjix);
end
