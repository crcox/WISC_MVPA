function summarize_jobs()
	allfiles = dir('./');
	dirs = allfiles([allfiles.isdir]);
	jobdirs = dirs(cellfun(@(x) ~(strcmp(x,'shared')||strcmp(x,'.')||strcmp(x,'..')), {dirs.name}));
	
	
end
