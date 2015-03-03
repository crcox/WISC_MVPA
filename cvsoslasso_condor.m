function [fitObj,fitObjRaw] = cvsoslasso_condor(X,Y,CV,LAMSET,ALPHA,GroupInfo,opts)
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

      % overlap_2stage returns matrices with one column per subject.
			[B(:,a:b),B_raw(:,a:b),W] = overlap_2stage(Ytrain,Xtrain,GroupInfo,LAMBDA,ALPHA,opts);
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
    % Data model:
    % Each subject is an element of a structured array.
    fitObj(nsubj,ncv) = struct('betas',[],'a0',[],'Y',[],'Yh',[],'lambda',[],'alpha',[],'omit',[]);
		SubjectLabels = mod(0:(size(B,2)-1),nsubj)+1;
		for ss = 1:nsubj
      for cc = 1:ncv
				svox = [false;GroupInfo.SubjectMask(:,ss)]; % adjusted for bias unit. 
        z = SubjectLabels == ss;
        fitObj(ss,cc).a0 = B(1,z);
        fitObj(ss,cc).betas = B(svox,z);
				fitObj(ss,cc).Y = Y{ss};
        fitObj(ss,cc).Yh = full(X{ss} * B(:,z));
        fitObj(ss,cc).lambda = LAMSET;
        fitObj(ss,cc).alpha = ALPHA;
        fitObj(ss,cc).testset = CV{ss}(:,cc);
      end
		end
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

