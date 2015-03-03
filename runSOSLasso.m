function [fitObj,GroupInfo,metadata, X] = main()
	% This contains the following:
	% regressors coding out the different conditions that I might want to
	% classify. I called these TrueFaces, TruePlaces, TrueThings, and
	% you'll see reference to those below when I define Y and the CV
	% blocks.
	%
	% A logical matrix coding out my CV blocks. There is a row for each
	% trial/TR/whatever, and a column for each CV block. For me, true means
	% "omit this TR on this run of the analysis."
	%
	% ***xyz_tlrc*** THIS IS IMPORTANT. This is a cell array, with one cell
	% per subject, that contains the tlrc warped coordinates for each
	% subject. Here, we are doing something a little unusual, so I want to
	% be clear: the functional data are NOT warped, per se. We just warped
	% the coordinates themselves. To obtain these weights we ran a warping
	% procedure on the anatomicals, and the extracted the transformation
	% matrix, and applied that manually to each subject's coordinates.
	% SOSLasso will form groups based on these coordinates. It involves
	% defining a space that is large enough to contain all the points for
	% all the subjects. It is possible that there are gaps between voxels
	% after this transformation, and that voxels don't exactly line up
	% across subjects. THAT IS OK. SOSLasso doesn't need them to be
	% perfectly aligned. Right now, there is a procedure built into
	% prep_data() (see below) that will generate groups based on spatial
	% proximity. You specify in the parameters the group size (either in
	% three dimensions, or provide a scalar for isotropic groups), and it
	% will basically divide to common space up into however many groups it
	% takes to panel the whole space. You also specify GroupShift, which
	% will say how much to offset each group---groups can overlap, so
	% GroupShift=1 means that you will overlap a lot and
	% GroupShift=GroupSize means that there is no overlap.
	%
	% Of course, feel free to add whatever else you want to metadata.mat.

	% This is the variable I am using to specify which part of the data to
	% omit entirely, for use as a final hold out set. Cross validation over
	% the remaining blocks of performed to find the best alpha and lambda.

	params = loadjson('params.json');
	load(params.metadata);
	nSubjects = length(metadata);
	for i = 1:nSubjects
		if i == 1
			ncv_total = size(metadata(i).CVBLOCKS,2);
		else
			prev = ncv_total;
			ncv_total = size(metadata(i).CVBLOCKS,2);
			if prev ~= ncv_total
				error('Subjects have different number of CVBLOCKS.')
			end
		end
	end

	if iscellstr(params.data)
		nDatafiles = length(params.data);
		if nDatafiles ~= nSubjects
			error('nDatafiles does not equal nSubjects.')
		end
		X = cell(nDatafiles, 1);
		for i = 1:length(params.data)
			temp = load(params.data{i},'X');
			X{i} = temp.X;
		end
		clear temp;
	else
		load(params.data, 'X')
		if ~iscell(X)
			X = cell(X);
		end
	end

	OMIT = params.FinalHoldoutSet;

	% If you do init_opts([]) you will get a template structure. Running
	% init_opts(struct) will fill in anything essential that you might have
	% missed, and perform some quality assurance.
	params = init_opts(params);

	%% Set up Y
	Y = {metadata.(params.TargetCategory)}';

	%% Project subjects' data in to Shared Space and create GroupInfo.
	% We do not want the intercept unit to belong to a group, so add it after
	% prepping all the group info.
	[X,GroupInfo] = prep_data(X,metadata,params);
	for ss = 1:nSubjects
		if params.UseIntercept == true
			X{ss} = [ones(size(X{ss},1),1),X{ss}]; % fit an intercept.
		else
			X{ss} = [zeros(size(X{ss},1),1),X{ss}]; % do not fit an intercept.
		end
	end

	ALPHA = params.alpha; % 0 lasso, 1 group lasso
	if isfield(params,'lambda') && ~isempty(params.lambda)
		LAMSET = params.lambda; % Scales overall regularization strength
	else
		LAMSET = lambdaseq(X,Y,ALPHA,100);
	end

	%% CV selections
	test = cellfun(@(x) x(:,OMIT), {metadata.CVBLOCKS}, 'unif', 0);
	Xtrain = cell(size(X));
	Ytrain = cell(size(Y));
	CV = cell(size(Y));
	for i = 1:nSubjects
		Xtrain{i} = X{i}(~test{i},:);
		Ytrain{i} = Y{i}(~test{i});
		CV{i} = metadata(i).CVBLOCKS(~test{i},(1:ncv_total)~=OMIT);
	end

 	[fitObj,fitObj_raw] = cvsoslasso_condor(Xtrain,Ytrain,CV,LAMSET,ALPHA,GroupInfo,params);
%		 save(sprintf('jlp+soslasso_mu%02d_cv_%02d.mat',WhichMu,OMIT),'fitObj','fitObjRaw','OMIT','WhichMu');
%		 [fitObj, fitObjRaw] = deal(0);

  [fitObj.omit] = deal(OMIT);
  save('fitObj.mat','fitObj');
% 	writeResults('results',fitObj)
% 	if params.SaveRawResults == true
%    [fitObj_raw.omit] = deal(OMIT);
% 		writeResults('results_NoDebiasing',fitObj_raw);
% 	end

	f=fopen('ALL_DONE','w');
	fclose(f);

	function writeResults(outdir,fitObj)
    % FitObj is subject x cv structured array

    % First, save the full fitObj
		mkdir(outdir);
		save(fullfile(outdir,'fitObj.mat'),'fitObj');

    % Then create a directory for each subject, 

		pathBin = @(x,y,z) fullfile(x,sprintf('%s_%02d.bin',y, z));
		save(fullfile(outdir,'fitObj.mat'),'fitObj');
		for i = 1:nSubjects
			writeBinMatrix(pathBin(outdir,'beta',i), full([fitObj.a0{i};fitObj.betas{i}]));
			writeBinMatrix(pathBin(outdir,'fittedVals',i), full(fitObj.Yh{i}));
			writeBinMatrix(pathBin(outdir,'dprime_test',i), full(fitObj.dp{i}));
			writeBinMatrix(pathBin(outdir,'dprime_train',i), full(fitObj.dpt{i}));
		end
		cvind = 1:ncv_total;
		cvind(OMIT) = [];
		[cv,lam,subj,train,omit,alpha] = ndgrid(cvind,LAMSET,1:nSubjects,[1,0],OMIT,ALPHA);
		dpvec = [cell2mat(cellfun(@(x) x(:), fitObj.dp,'unif',false)); ...
						 cell2mat(cellfun(@(x) x(:), fitObj.dpt,'unif',false))];
		dptbl = [cv(:),lam(:),subj(:),train(:),omit(:),alpha(:),dpvec];
		dptbl_header = {'cv','lam','subj','test','omit'};
		save(fullfile(outdir,'dprime_table.mat'),'dptbl','dptbl_header');
		csvwrite(fullfile(outdir,'dprime_table.csv'),dptbl);
		writeBinTable(fullfile(outdir,'dprime_table.bin'),dptbl,'compress',false);
		% NB Cannot compress on CONDOR since that would require java.
	end
end

function writeBinMatrix(filename,X)
% Column-major order.
	f = fopen(fullfile(filename),'w');
	fwrite(f,size(X), 'int', 'ieee-le'); % 4 byte
	fwrite(f, X, 'double', 'ieee-le');   % 8 byte
	fclose(f);
end
function writeBinTable(filename,tbl,varargin)
% Column-major order
% Cannot compress on CONDOR since that would require java.
	p = inputParser;
	p.addOptional('compress',false,@islogical);
	parse(p, varargin{:});
	f = fopen(filename,'w');
	fwrite(f,size(tbl), 'uint', 'ieee-le');    % 4 byte, dims
	fwrite(f, tbl(:,1), 'uint8', 'ieee-le');   % 1 byte, cv
	fwrite(f, tbl(:,2), 'float', 'ieee-le');   % 4 byte, lam
	fwrite(f, tbl(:,3), 'uint8', 'ieee-le');   % 1 byte, subj
	fwrite(f, tbl(:,4), 'uint8', 'ieee-le');   % 1 byte, isTest
	fwrite(f, tbl(:,5), 'uint8', 'ieee-le');   % 1 byte, omit
	fwrite(f, tbl(:,6), 'float', 'ieee-le');   % 4 byte, alpha
	fwrite(f, tbl(:,7), 'float', 'ieee-le');   % 4 byte, dprime
	fclose(f);
	if p.Results.compress == true
		gzip(filename);
		delete(filename)
	end
end
