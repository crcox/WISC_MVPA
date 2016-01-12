function WriteVoxels(outdir,results,metadata,varargin)
% This function will write out a single directory worth of data, where a
% directory will include all cross validation instances for all subjects for a
% particular configuration of parameters.
% The file naming convention is subject_finalholdout_cvholdout.CoordLabel.
  p = inputParser();
  addRequired(p, 'outdir');
  addRequired(p, 'results');
  addRequired(p, 'metadata');
  addParameter(p, 'CoordLabel', 'mni');
  addParameter(p, 'DataVar','Wz');
  addParameter(p, 'Method','soslasso');
  addParameter(p, 'SubjectNumberFMT', 's%d_rep.mat');
  parse(p, outdir, results, metadata, varargin{:});

  outdir    = p.Results.outdir;
  results   = p.Results.results;
  xyzlab    = p.Results.CoordLabel;
  datavar   = p.Results.DataVar;
  method    = p.Results.Method;

  n = length(results);
  fprintf('Writing files:\n');
  for i = 1:n
    R  = results(i);
    sj = R.subject;
    fh = R.finalholdout;
    cv = R.cvholdout;
    M = metadata(sj);
    idx = find(strcmp(xyzlab,{M.coords.orientation}));
    z = M.filter(2).filter;
    xyz = M.coords(idx).xyz(z,:);
    if isfield(R,'Wix')
      vx = R.Wix;
    else
      vx = (1:size(xyz,1))';
    end


    if length(cv) > 1
      fname = sprintf('%02d_%02d_xx.%s', sj, fh, xyzlab);
    else
      fname = sprintf('%02d_%02d_%02d.%s', sj, fh, cv, xyzlab);
    end

    if strcmp(method, 'searchlight');
      ra = R.slradius;
      fname = sprintf('%02d_%02d_%02d.%s', sj, fh, ra, xyzlab);
    end
    fpath = fullfile(outdir, fname);
    fprintf('%s\n', fpath);

    wz  = R.(datavar)(:);
    dlmwrite(fpath, [xyz(vx,:),wz], ' ');
  end
  fprintf('Done.\n');
end
