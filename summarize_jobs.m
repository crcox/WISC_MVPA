function S=summarize_jobs()
	jdat = loadjson('master.json');
	nlam = length(jdat.lambda);
	
	%% Arrange the free parameters
	% ASSUMPTION: A common set of variable are being varied over configs.
	% This means it is possible to collapse the cell array of structures read
	% in from the master.json into a 1 x n structured array for convenience.
	for i=1:length(jdat.config);
		config(i) = jdat.config{i}; %#ok<AGROW>
	end
	
	% ASSUMPTION: The variables that can be varied are:
	% - lambda (but this is done within jobs, currently.
	% - alpha
	% - GroupShift
	% - GroupSize
	PLabs = {'alpha','GroupShift','GroupSize'};
	P = zeros(length(config),3);
	sz = zeros(1,3);
	if isfield(config,'alpha')
		P(:,1) = [config.alpha];
		sz(1) = length(unique([config.alpha]));
	end
	if isfield(config,'GroupShift')
		P(:,2) = [config.GroupShift];
		sz(2) = length(unique([config.GroupShift]));
	end
	if isfield(config,'GroupSize')
		P(:,3) = [config.GroupSize];
		sz(3) = length(unique([config.GroupSize]));
	end
	PLabs = {'lambda',PLabs(any(P))};
	sz = [length(jdat.lambda),sz(any(P))];
	P = P(:,any(P));
	
	% LOGIC: This sortrows() call is going to return indexes that will allow
	% looping over the configs in a predictable way with respect to alpha,
	% GroupShift, and GroupSize.  
	[~,ix] = sortrows(P,size(P,2):-1:1);
	
	S = zeros(sz);
	for i = 1:length(config)
		% NOTE: ix(i) is ensuring that we load job results in the right order
		% so the summary stats can be predicably stored.
		jobdir = sprintf('%03d',ix(i)-1);
		load(fullfile(jobdir,'fitObj.mat'));
		[nsubj, ncv] = size(fitObj);
		dps = zeros(nsubj, nlam);
		for ss = 1:nsubj
			tmp = zeros(ncv, nlam);
			for cc = 1:ncv
				tmp(cc,:) = dprimeCV(fitObj(ss,cc).Y, fitObj(ss,cc).Yh, fitObj(ss,cc).testset);
			end
			dps(ss,:) = mean(tmp);
		end
		b = nlam * i;
		a = b - (nlam-1);
		S(a:b) = mean(dps);
	end
	
	[~,maxind] = max(S(:));

	% NOTE: Note the use of ix(...). This maps an index into S to a job
	% index.
	% LOGIC: The first dimension is lambda, which is the same for each job.
	% By integer-dividing by nlam (and rounding up), we essentially get a
	% column index into reshape(S,sz(1),prod(sz(2:end)). Those column indexes
	% map to individual jobs, and ix() contains the mapping to go from one of
	% those indexes to a job index.
	MaxCfgInd = ix(idivide(maxind,uint32(nlam),'ceil'));
	MaxLamInd = rem(maxind,uint32(nlam));
	
	fprintf('Best parameters (Job %d):\n', MaxCfgInd-1);
	fprintf('\talpha: %.3f\n', config(MaxCfgInd).alpha);
	fprintf('\tlambda: %.3f\n', jdat.lambda(MaxLamInd));
	fprintf('\tGroupSize: %d\n', config(MaxCfgInd).GroupSize);
	fprintf('\tGroupShift: %d\n', config(MaxCfgInd).GroupShift);
	fprintf('dprime: %.3f\n', S(maxind));
	
	BestCfg = rmfield(jdat,'config');
	BestCfg.lambda = jdat.lambda(MaxLamInd);
	if isfield(config,'alpha')
		BestCfg.alpha = config(MaxCfgInd).alpha;
	end
	if isfield(config,'GroupShift')
		BestCfg.GroupShift = config(MaxCfgInd).GroupShift;
	end
	if isfield(config,'GroupSize')
		BestCfg.GroupSize = config(MaxCfgInd).GroupSize;
	end
	
	disp(BestCfg);
	savejson('',BestCfg,'bestcfg.json');
end
