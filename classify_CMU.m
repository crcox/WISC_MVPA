clear
clc;
close all;

% load the data
load CMU_FMRI_551
doglasso = false;


g = G{end};
if isempty(g)
    G = G(1:end-1);
end
fprintf('FMRI DATA LOADED \n');
X1 = Xsure{1};
numpersons = length(Xsure);
numsamples = size(X1,1);
numvoxels = size(X1,2);
clear X1;

% set if doing non overlapping group lasso
if doglasso
   % define group lasso groups
   G = cell(numvoxels,1);
   for ii = 1:numvoxels
      G(ii) = {ii};  
   end
end


X = cell(0);
% normalize CMU data
for person = 1:numpersons
    Xtemp = Xsure{person};
    Xtemp = Xtemp - repmat(min(Xtemp),size(Xtemp,1),1); % mean center
    s = std(Xtemp);
    iszer = find(s==0); s(iszer) = 1; % will not be dividing by 0
    X = [X; {[Xtemp]}];
end
clear Xsure;
fprintf('Data normalized. Mean centered and unit std. dev.  \n')


%% problem parameters

% test set = 8*32
% training = 24*32
% validating = 8*32
% total = 1280

lamset = [0.001, 0.01, 0.1, 0.5, 1, 5, 10]; %linspace(1e-4,1e-1,25); % regularizer values
gamset = 2.^[-3:3];

% shuffle only trials
shuffletrial = randsample(40,40);
shuffle = [];
for ii = 1:40
    i = shuffletrial(ii);
    shuffle = [shuffle [1+(i-1)*32:i*32]];
end

testinds = shuffle(1:8*32);
stest = randsample(length(testinds),length(testinds));
testinds = testinds(stest);
numcvs = 4; cvsize = 8*32;
Beta_CV = cell(0);
tic
fprintf('beginning cv ... \n')

%% cross validation module
for cv = 1:numcvs
    
    % make training and hold out X and Y;
    cvind = (cv-1)*cvsize+1:(cv)*cvsize;
    cvind = cvind + length(testinds);
    trainind = setdiff(shuffle(length(testinds)+1:numsamples),shuffle(cvind));
    scv = randsample(length(cvind),length(cvind));
    cvind = cvind(scv);
    str = randsample(length(trainind),length(trainind));
    trainind = trainind(str);
    trainX = cell(0);trainY = cell(0);testX = cell(0); testY = cell(0);
    for person = 1:numpersons
        Xtemp = X{person};
        Ytemp = Y{person};
        Xtempt = Xtemp(trainind,:); %training data
        Xtempc = Xtemp(shuffle(cvind),:);    % CV data
        Ytempt = Ytemp(trainind,:);
        Ytempc = Ytemp(shuffle(cvind),:);
        trainX = [trainX ; {[Xtempt]}];
        testX  = [testX ; {[Xtempc]}];
        trainY = [trainY ; {[Ytempt]}];
        testY  = [testY ; {[Ytempc]}];
    end
    
    % run the method over all lambda values
    [Xo, groups, group_arr] = makeA_multitask(trainX,G);
    
    fprintf('matrices replicated \n')
%     % add in a bias term for each subject
%     for person = 1:numpersons
%         
%         Xtemp = Xo{person};
%         Xtemp = [ones(size(Xtemp,1),1) Xtemp];
%         Xo{person}= Xtemp;
%         
%     end
%     fprintf('Bias Term Added \n')
    
    y_for_err = repmat(Ytempc,1,numpersons);
    a = y_for_err(:);
    lamind = 0;
    for lam = lamset
        gamind = 0;
        lamind = 1+lamind;
        
        for gam = gamset
            gamind = 1+gamind;
            
            
            [Betahat,~] = overlap_2stage_CMU(trainY,Xo,X,G,group_arr,groups, lam,gam);
            Betahat = sparse(Betahat);
            
            %cross validate

            clear b;
            for person = 1:numpersons
                B = Betahat(:,person);
                x_for_err = testX{person};
                x_for_err = [ones(size(x_for_err,1),1)  x_for_err];
                b(:,person) = sign(x_for_err*B);
            end
            b = b(:);
            L = length(a);
            err = sum(sign(a)~=sign(b))/L;
            
            disp(err);
            disp(nnz(Betahat)/numel(Betahat));
            err = full(err);
            CV_ERROR(cv,lamind, gamind) = err;
            
            fprintf('.')
        end
        
        toc
        fprintf(2,'%d out of %d Lambda done. CV = %d \n',lamind, length(lamset), cv);
    end
    save CMU_DIAG CV_ERROR lamind gamind
end
fprintf('\n ALL CV DONE \n')
ALLERRS = CV_ERROR;
CV_ERROR = mean(CV_ERROR,1);

%% PICK THE BEST REGULARIZATION PARAMETER AND RELEARN MODEL
CV_ERROR = reshape(CV_ERROR,[size(CV_ERROR,2),size(CV_ERROR,3)]);
[row,col] = find(CV_ERROR == min(min(CV_ERROR)));

fprintf('min. training error = %f \',min(min(CV_ERROR)));

row = row(end); col= col(end);
lam = lamset(row);  % pick the lambda that minimzes the CV error
gam = gamset(col);
trainind = shuffle(1+length(testinds):numsamples);  %entire training set
trainX = cell(0);trainY = cell(0);
for person = 1:numpersons
    Xtemp = X{person};
    Ytemp = Y{person};
    Xtempt = Xtemp(trainind,:);
    Ytempt = Ytemp(trainind,:);
    trainX = [trainX ; {[Xtempt]}];
    trainY = [trainY ; {[Ytempt]}];
end

fprintf('replicating all data again \n')
[Xo, groups, group_arr] = makeA_multitask(trainX,G);


fprintf('training final model ...')
[Betahat,~] = overlap_2stage_CMU(trainY,Xo,X,G,group_arr,groups, lam,gam);
        
Betahat  =sparse(Betahat);
fprintf('done \n')
% Betahat is the final model

%% TEST MODEL PERFORMANCE ON TEST SET

finX = cell(0); finY = cell(0);
for person = 1:numpersons
    Xtemp = X{person};
    Ytemp = Y{person};
    Xtemp = Xtemp(testinds,:);
    Ytemp = Ytemp(testinds,:);
    finX = [finX ; {[Xtemp]}];
    finY = [finY ; {[Ytemp]}];
end
y_for_err = repmat(Ytemp,1,numpersons);
a = y_for_err(:);
clear b;
for person = 1:numpersons
    B = Betahat(:,person);
    x_for_err = finX{person};
    x_for_err = [ones(size(x_for_err,1),1)  x_for_err];

    b(:,person) = sign(x_for_err*B);
end
b = b(:);
L = length(a);
fprintf('overall error = ')
err = sum(sign(a) ~= sign(b))/L %this is the overall error
%%
c = clock;
c = c(2:end-1); % month day hour min
strin = strcat(num2str(c(1)),num2str(c(2)),'_',num2str(c(3)),num2str(c(4)));

str = strcat('save CMU_SOSLASSO_',strin,' Betahat X Y shuffle err ALLERRS lam gam');
eval(str);
        
fprintf('\n DATA SAVED \n')