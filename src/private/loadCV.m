function cvind = loadCV(cvpath, cv_var, cvscheme, rowfilter)
  % In this dataset, trials are ordered the same way across subjects, so each
  % subject gets a copy of the same CV scheme.
  rowfilter = ascell(rowfilter);
  N = numel(rowfilter);
  tmp = load(cvpath, cv_var);
  CV = tmp.(cv_var); clear tmp;
  if N > 1
    tmp = {CV(:, cvscheme)};
    cvind = repmat(tmp, N, 1);
    for i=1:N
      cvind{i} = cvind{i}(rowfilter{i});
    end
  else
    rowfilter = uncell(rowfilter);
    cvind = CV(rowfilter, cvscheme);
  end
end
