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
  addParameter(p, 'SubjectNumberFMT', 's%d_rep.mat');
  parse(p, outdir, results, metadata, varargin{:});

  outdir    = p.Results.outdir;
  results   = p.Results.results;
  xyzlab    = p.Results.CoordLabel;

  n = length(results);
  fprintf('Writing files:\n');
  for i = 1:n
    R  = results(i);
    sj = R.subject;
    fh = R.finalholdout;
    cv = R.cvholdout;
    vx = R.Wix;

    M = metadata(sj);

    fname = sprintf('%02d_%02d_%02d.%s', sj, fh, cv, xyzlab);
    fpath = fullfile(outdir, fname);
    fprintf('%s\n', fpath);

    idx = find(strcmp(xyzlab,{M.coords.orientation}));
    xyz = M.coords(idx).xyz;
    wz  = R.Wz(:);
    dlmwrite(fpath, [xyz(vx,:),wz], ' ');
  end
  fprintf('Done.\n');
end
