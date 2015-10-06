function filepath = updateFilePath(file, path, suffix)
  file     = ascell(file);
  N        = length(file);
  filepath = cell(1,N);
  for i = 1:N
    [~,f,e] = fileparts(file{i});
    if nargin>2
      f = sprintf('%s_%s%s',f,suffix,e);
    else
      f = sprintf('%s%s',f,e);
    end
    filepath{i} = fullfile(path,f);
  end
  if N == 1
    filepath = filepath{1};
  end
end
