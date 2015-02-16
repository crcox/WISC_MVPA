function [fitObj,fitObjRaw] = cvsoslasso_condor(X,Y,CV,LAMSET,MU,GroupInfo,opts)
% X : Nx1 Cell array, where N is the number of subjects.
% Y : Nx1 Cell array, where N is the number of subjects.

	nsubj = length(X);
	nlam = length(LAMSET);
	ncv = size(CV{1},2);

	if exist(fullfile(pwd, 'CHECKPOINT.mat'), 'file') == 2
		[B,B_raw,cv_start,lam_start,opts] = loadFromCheckpoint();
	else
		opts = init_opts(opts);
		cv_start = 1;
		lam_start = 1;
		[B,B_raw] = deal(zeros(size(X{1},2),ncv*nlam*nsubj));
	end

	for cc = cv_start:ncv
		fprintf('cv: %d, lam:', cc);
		[Xtrain, Ytrain] = selectTrainingSets();

		for ll = lam_start:nlam
			fprintf(' %d', ll);
			LAMBDA = LAMSET(ll);
			%% Fit soslasso model to the training set.
			a = sub2ind([nsubj,ncv,nlam],1,cc,ll);
			b = sub2ind([nsubj,ncv,nlam],nsubj,cc,ll);

			[B(:,a:b),B_raw(:,a:b),W] = overlap_2stage(Ytrain,Xtrain,GroupInfo,LAMBDA,MU,opts);
			opts.W0=W(2:end,:); % Warm start
			opts.C0=W(1,:); % Warm start
			opts.init=2; % set to 1 to let warm start happen.
			B = sparse(B);
			B_raw = sparse(B_raw);
			save('CHECKPOINT.mat','B','B_raw','cc','ll','opts');
		end
		fprintf('\n');
		opts = rmfield(opts,'W0');
		opts = rmfield(opts,'C0');
		opts.init=2;
		lam_start = 1;
	end

	%% Evaluate and store debiased solutions
	fitObj = evaluateAndStoreModel(B);

	%% Evaluate and store raw solutions
	fitObjRaw = evaluateAndStoreModel(B_raw);

	delete(fullfile(pwd,'CHECKPOINT.mat'));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%             BEGIN NESTED FUNCTIONS                 %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	function [Xtrain, Ytrain] = selectTrainingSets()
		[Xtrain, Ytrain] = deal(cell(nsubj,1));
		for ss = 1:nsubj
			Xtrain{ss} = X{ss}(~CV{ss}(:,cc),:);
			y = double(Y{ss}(~CV{ss}(:,cc)));
			y(y==0) = -1;
			Ytrain{ss} = y;
		end
	end

	function fitObj = evaluateAndStoreModel(B)
		betas = cell(size(X));
		a0 = cell(size(X));
		Yh = cell(size(X));
		dp = cell(size(X));
		dpt = cell(size(X));
		SubjectLabels = mod(0:(b-1),nsubj)+1;
		for ss = 1:nsubj
			z = SubjectLabels == ss;
			a0{ss} = B(1,z);
			betas{ss} = B(2:end,z);
			Yh{ss} = full(X{ss} * B(:,z));
			dp{ss} = reshape(dprimeCV(Y{ss}>0,Yh{ss}>0, CV{ss}), ncv, nlam);
			dpt{ss} = reshape(dprimeCV(Y{ss},Yh{ss}, ~CV{ss}), ncv, nlam);
		end

		[~,bestLambdaInd] = max(mean(cell2mat(dp)));
		bestLambda = LAMSET(bestLambdaInd);

		fitObj.betas = betas;
		fitObj.a0 = a0;
		fitObj.lambda = LAMSET;
		fitObj.mu = MU;
		fitObj.bestLambda = bestLambda;
		fitObj.bestLambdaInd = bestLambdaInd;
		fitObj.dp = dp;
		fitObj.dpt = dpt;
		fitObj.Yh = Yh;
		fitObj.Y = Y;
		fitObj.CV = CV;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             BEGIN PRIVATE FUNCTIONS                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B,B_raw,cc,ll,opts] = loadFromCheckpoint() %#ok<STOUT>
%% Fit soslasso model to the training set.
	load('CHECKPOINT.mat');
	fprintf('Resuming from cv: %d, lam: %d.\n',cc,ll);
	% 'B','B_raw','cc','ll','opts'
end

