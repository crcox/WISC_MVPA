function [dp,counts] = dprime(truth,prediction,jitter)
% Takes two binary vectors as input, one this represents the "truth", and
% the other the prediction of some classifier.
if nargin == 2
    jitter = 0.0005;
end
T = truth > 0;
P = prediction > 0;

hit = sum(bsxfun(@and,T,P));
fa = sum(bsxfun(@and,~T,P));

hitr = hit./sum(T);
far = fa./sum(~T);
if isnan(hitr)
    error('There are no targets in the test set.')
elseif isnan(far)
    error('The test set is entirely targets.')
end

hitr(hitr==0) = jitter;   far(far==0) = jitter;
hitr(hitr==1) = 1-jitter; far(far==1) = 1-jitter;

dp = norminv(hitr) - norminv(far);

hitr(hitr==jitter) = 0;   far(far==jitter) = 0;
hitr(hitr==(1-jitter)) = 1; far(far==(1-jitter)) = 1;

counts.targets = sum(T);
counts.distractors = sum(~T);
counts.hits = hit;
counts.hitrate = hitr;
counts.falsealarms = fa;
counts.falsealarmrate = far;
