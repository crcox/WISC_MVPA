function id = extractSubjectID(datafile)
  [~,fname,~] = fileparts(datafile);
  str = regexp(fname,'[0-9]+','match');
  id = sscanf(str{1},'%d');
  if isempty(id)
    error('Failed to extract subject id from data filename %s. Exiting...', datafile);
  end
end
