%
% compute autocorrelation inside each column, at various lags
%
% in:
% - data matrix - #observations x #variables
% - lags
%
% out:
% - autocorrelation - #lags x #variables
%
% notes:
% - ideally, mean and var of each feature would be computed once for the smallest lag and then
% updated for each lag by removing just the unwanted values; for now, it's simpler to just recompute
%
% history:
% 20100321 - created - fpereira@princeton.edu
%

function [autocorrelation] = computeAutocorrelation( varargin )

if nargin < 2
  fprintf('syntax: computeAutocorrelation(matrix,lags)\n');return;
end

matrix = varargin{1}; [n,m] = size(matrix);
lags   = varargin{2}; nl = length(lags);

if sum(lags < 0); fprintf('error: lags must be >= 0\n'); return; end
if sum(lags > (n-2)); fprintf('error: lags must be <= #observations-2\n');return;end

autocorrelation = zeros(nl,m);

for ip = 1:nl
  lag = lags(ip);
  
  rangeA = 1:(n-lag);
  rangeB = (1+lag):n;
  ncut = n-lag;

  % crop the matrices (B is original shifted by lag, A is original - lag at the end)
  tmpA = matrix(rangeA,:);
  tmpB = matrix(rangeB,:);
  
  % normalize them so that each column is mean 0 stdv 1
  meanA = mean(tmpA,1); stdvA = std( tmpA,0,1);
  meanB = mean(tmpB,1); stdvB = std( tmpB,0,1);
  
  tmpA = tmpA  - repmat(meanA,ncut,1); tmpA = tmpA ./ repmat(stdvA,ncut,1);
  tmpB = tmpB  - repmat(meanB,ncut,1); tmpB = tmpB ./ repmat(stdvB,ncut,1);

  % and compute the correlation between each column
  autocorrelation(ip,:) = sum(tmpA .* tmpB,1) / (ncut-1);
end

function [] = testThis()

% each example depends on previous one
n = 100; m = 100;
X = randn(n,m);
X(:,round(m/2):m) = 0;
Xcorrelated = X;

a = [0.25 0.25 0.5];
la = length(a);

for ip = la:n
  range = (ip-la+1):ip;
  Xcorrelated(ip,:) = sum(repmat(a',1,m) .* X(range,:),1);
end
Xcorrelated = Xcorrelated + randn(size(Xcorrelated))*0.5;
lags = [1 2 3 4 5 6];
[autocorrelation] = computeAutocorrelation( Xcorrelated, lags );

tmp = repmat(n - lags',1,m);

tauto = autocorrelation ./ sqrt((1-autocorrelation.^2)./(tmp-2));
tsig  = zeros(size(tauto));
for ip = 1:length(lags); tsig(ip,:) = 1-tcdf(tauto(ip,:),tmp(ip,1)-2); end
[discard,threshold] = computeFDR(tsig(:),0.01);
clf;imagesc(tsig<=threshold);pause

zauto = (log((1+autocorrelation)./(1-autocorrelation))/2) .* sqrt(tmp-1);
zsig = zeros(size(zauto));
for ip = 1:length(lags); zsig(ip,:) = 1-normcdf(zauto(ip,:)); end
[discard,threshold] = computeFDR(zsig(:),0.01);
clf;imagesc(zsig<=threshold);pause

clf;
subplot(3,1,1);hist(autocorrelation(1,:),20);axis([-1 1 0 Inf]);
subplot(3,1,2);hist(autocorrelation(2,:),20);axis([-1 1 0 Inf]);
subplot(3,1,3);hist(autocorrelation(3,:),20);axis([-1 1 0 Inf]);



% forward, BOLD style
x = 0:12; bold = 3*gampdf(x-2,3,1);
ntrials = 50;
trialsequences{1} = [1 1 1 0 0 0 0 0 0 0 0]
m = 20; nblank = 0;
ns = length(trialsequences);
boldlen = length(bold);

for it = 1:ns
  trialsequence = trialsequences{it};
  duration = length(trialsequence);
  
  n = duration * ntrials;
  nplus = n + boldlen - 1;

  timecourse = repmat(trialsequence',ntrials,1); 
  X = repmat(timecourse,1,m);
  convolved = zeros(nplus,1);
  for i = 1:n
    if timecourse(i)
      range = i:(i+boldlen-1);
      convolved(range) = convolved(range) + bold';
    end
  end
  Xconvolved = repmat(convolved,1,m) + randn(nplus,m)*0.5;

  lags = 1:15;
  [autocorrelation] = computeAutocorrelation( Xconvolved, lags );

  imagesc(autocorrelation,[-1 1]);colorbar('vert');
  pause
end




Xcorrelated = zeros(n,m);
a = [1 2 3 4 5 6 4 2]; a = a/sum(a);

la = length(a);

for ip = 1:(n-la)
  range = (ip+1):(ip+1+la-1);
  Xcorrelated(ip,:) = Xcorrelated(ip,:) + sum(repmat(a',1,m) .* X(range,:),1);
end

clf;subplot(1,2,1);imagesc(X);subplot(1,2,2);imagesc(Xcorrelated);

% with real data

study = 'data-twocategories';
subject = '01481B';
%preprocessingSteps{1} = 'voxelRunMeanZero';
%preprocessingSteps{2} = 'voxelRunStdvOne';
preprocessingSteps = {'imageMeanZeroStdvOne'};
[meta,examples,C,studyInfo] = datasetLoader(study,subject,{},preprocessingSteps);
lags = [1 2 3 4 5 6 7 8 9 10];
[autocorrelation] = computeAutocorrelation( examples, lags );

examplesNoise = randn(size(examples));
[autocorrelationNoise] = computeAutocorrelation( examplesNoise, lags );

clf;subplot(2,1,1);hist(autocorrelation(1,:),20);axis([-1 1 0 Inf]);
subplot(2,1,2);hist(autocorrelation(2,:),20);axis([-1 1 0 Inf]);
subplot(2,1,2);hist(autocorrelationNoise(1,:),20);axis([-1 1 0 Inf]);

ranksum(autocorrelation(1,:),autocorrelationNoise(1,:))
