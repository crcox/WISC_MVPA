function [S,NZ]=summarize_jobs()
	jdat = loadjson('master.json');
	%% Arrange the free parameters
	% ASSUMPTION: A common set of variable are being varied over configs.
	% This means it is possible to collapse the cell array of structures read
	% in from the master.json into a 1 x n structured array for convenience.
	for i=1:length(jdat.config);
		config(i) = loadjson(fullfile(sprintf('%03d',i-1),'params.json')); %#ok<AGROW>
	end

	% ASSUMPTION: The variables that can be varied are:
	% - lambda
	% - alpha
	% - GroupShift
	% - GroupSize
	PLabs = {'lambda','alpha','GroupShift','GroupSize'};
	P = zeros(length(config),length(PLabs));
	sz = zeros(1,length(PLabs));
	if isfield(config,'lambda')
		P(:,1) = [config.lambda];
		lamset = unique([config.lambda]);
		sz(1) = length(lamset);
		nlam = length(lamset);
	end
	if isfield(config,'alpha')
		P(:,2) = [config.alpha];
		alphaset = unique([config.alpha]);
		sz(2) = length(alphaset);
	end
	if isfield(config,'GroupShift')
		P(:,3) = [config.GroupShift];
		shiftset = unique([config.GroupShift]);
		sz(3) = length(shiftset);
	end
	if isfield(config,'GroupSize')
		P(:,4) = [config.GroupSize];
		sizeset = unique([config.GroupSize]);
		sz(4) = length(sizeset);
	end
	PLabs = PLabs(any(P));
	P = P(:,any(P));

	% LOGIC: This sortrows() call is going to return indexes that will allow
	% looping over the configs in a predictable way with respect to alpha,
	% GroupShift, and GroupSize.
	[~,ix] = sortrows(P,size(P,2):-1:1);

	disp(sz);
	S = zeros(sz);
	NZ = zeros(sz);
	for i = 1:length(config)
		% NOTE: ix(i) is ensuring that we load job results in the right order
		% so the summary stats can be predicably stored.
		jobdir = sprintf('%03d',ix(i)-1);
		resultfile = fullfile(jobdir,'fitObj.mat');
		if exist(resultfile,'file');
			load(resultfile);
			fprintf('%s\n', resultfile)
		else
			fprintf('%s does not exist.\n', resultfile);
			continue
		end
		[nsubj, ncv] = size(fitObj);
		dps = zeros(nsubj,1);
		nzs = zeros(nsubj,1);
		for ll = 1:nlam
			for ss = 1:nsubj
				tmp = zeros(1,ncv);
				tmpnz = zeros(1,ncv);
				for cc = 1:ncv
					tmp(cc)   = dprimeCV(fitObj(ss,cc).Y, fitObj(ss,cc).Yh, fitObj(ss,cc).testset);
					tmpnz(cc) = sum(fitObj(ss,cc).betas~=0);
				end
				dps(ss) = mean(tmp);
				nzs(ss) = mean(tmpnz);
			end
		end
		S(i) = mean(dps);
		NZ(i) = mean(nzs);
	end

	[~,maxind] = max(S(:));

	% NOTE: Note the use of ix(...). This maps an index into S to a job
	% index.
	% LOGIC: The first dimension is lambda, which is the same for each job.
	% By integer-dividing by nlam (and rounding up), we essentially get a
	% column index into reshape(S,sz(1),prod(sz(2:end)). Those column indexes
	% map to individual jobs, and ix() contains the mapping to go from one of
	% those indexes to a job index.
% 	MaxCfgInd = ix(idivide(maxind,uint32(nlam),'ceil'));
% 	MaxLamInd = rem(maxind,uint32(nlam));
	MaxCfgInd = ix(maxind);

	fprintf('Best parameters (Job %d):\n', MaxCfgInd-1);
	fprintf('\talpha: %.3f\n', config(MaxCfgInd).alpha);
	fprintf('\tlambda: %.3f\n', config(MaxCfgInd).lambda);
	fprintf('\tGroupSize: %d\n', config(MaxCfgInd).GroupSize);
	fprintf('\tGroupShift: %d\n', config(MaxCfgInd).GroupShift);
	fprintf('dprime: %.3f\n', S(maxind));

	BestCfg = rmfield(jdat,'config');
	BestCfg.lambda = config(MaxCfgInd).lambda;
	if isfield(config,'alpha')
		BestCfg.alpha = config(MaxCfgInd).alpha;
	end
	if isfield(config,'GroupShift')
		BestCfg.GroupShift = config(MaxCfgInd).GroupShift;
	end
	if isfield(config,'GroupSize')
		BestCfg.GroupSize = config(MaxCfgInd).GroupSize;
	end
  BestCfg.isFinal = true;

	disp(BestCfg);
	savejson('',BestCfg,'bestcfg.json');
end
