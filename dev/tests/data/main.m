allfiles = dir('./');
dirs = allfiles([allfiles.isdir]);
jobdirs = dirs(cellfun(@(x) ~isempty(regexp(x,'[0-9]{3}')), {dirs.name}));

rootdir = pwd();
for i = 1:length(jobdirs);
	cd(jobdirs(i).name);
	allfiles=dir('./');
	files = allfiles(~[allfiles.isdir]);
	if ~any(strcmp('fitObj.mat',{files.name}));
		runSOSLasso();
	end
	cd(rootdir);
end