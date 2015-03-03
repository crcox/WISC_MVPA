allfiles = dir('./');
dirs = allfiles([allfiles.isdir]);
jobdirs = dirs(cellfun(@(x) ~isempty(regexp(x,'[0-9]{2}')), {dirs.name}));

rootdir = pwd();
for i = 1:length(jobdirs);
	cd(jobdirs(i).name);
	runSOSLasso();
	cd(rootdir);
end