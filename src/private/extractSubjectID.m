function id = extractSubjectID(datafile, FMT_subjid)
  [~,fname,~] = fileparts(datafile);
  id = sscanf(fname, FMT_subjid);
  if isempty(id)
    error('Failed to extract subject id from data filename %s. Exiting...', datafile);
  end
end
