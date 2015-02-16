%% Setup metadata for demo:
nSubjects = 10;
metadata(10).TrueFaces = [];
metadata(10).xyz_tlrc = [];
metadata(10).CVBLOCKS = [];
for ss = 1:nSubjects
	disp(ss)
	metadata(ss).CVBLOCKS = bsxfun(@eq,reshape(mod(0:999,10),100,10),0:9);
		metadata(ss).TrueFaces = [true(50,1);false(50,1)];
	while true
			metadata(ss).xyz_tlrc = unique([randsample(8,500,true),randsample(8,500,true),randsample(8,500,true)],'rows');
		if length(metadata(ss).xyz_tlrc)>250
			metadata(ss).xyz_tlrc = sortrows(metadata(ss).xyz_tlrc(randperm(length(metadata(ss).xyz_tlrc),250),:));
			break
		end
	end
end
save('metadata.mat','metadata');

%% Set up X
% Here is where you would load your data.
X = cell(nSubjects,1);
actVox = false(250,10);
for ss = 1:nSubjects
	s.X = randn(100,250);
    s.Y = metadata(ss).TrueFaces;
	metadata(ss).actVox = pdist2([4.5,4.5,4.5],metadata(ss).xyz_tlrc)<2;
	s.X(s.Y,metadata(ss).actVox) = s.X(s.Y,metadata(ss).actVox) + 3;
	X{ss} = bsxfun(@minus,s.X,mean(s.X));
	X{ss} = bsxfun(@rdivide,X{ss},std(X{ss}));
end
clear s;
save('alldata.mat', 'X');