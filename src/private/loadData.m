function [X,subjix] = loadData(datafile,data_var,rowfilter,colfilter,metadata, FMT_subjid)
  % Load data for multiple subjects, and apply filters.
  datafile  = ascell(datafile);
  rowfilter = ascell(rowfilter);
  colfilter = ascell(colfilter);
  N         = length(datafile);
  subjix    = zeros(1,N);
  X         = cell(N,1);
  for i = 1:N
    fprintf('Loading %s from  %s...\n', data_var, datafile{i});
    subjid    = extractSubjectID(datafile{i}, FMT_subjid);
    if ischar(subjid)
      subjix(i) = find(strcmp({metadata.subject}, subjid));
    else
      subjix(i) = find([metadata.subject] == subjid);
    end
    tmp       = load(datafile{i}, data_var);
    X{i}      = tmp.(data_var); clear tmp;
    X{i}      = X{i}(rowfilter{subjix(i)},colfilter{subjix(i)});
  end
end
