function [dp,err,counts] = dprimeCV(truth,prediction,cvblocks,jitter)
% Takes two binary vectors as input, one this represents the "truth", and
% the other the prediction of some classifier.
if nargin == 3
    jitter = 0.0005;
end

ntc = size(truth,2);
npc = size(prediction,2);
ncc = size(cvblocks,2);
[mm,mi] = max([ntc,npc,ncc]);
truth = repmat(truth,1,mm/ntc);
prediction = repmat(prediction,1,mm/npc);
cvblocks = repmat(cvblocks,1,mm/ncc);

T = truth > 0;
P = prediction > 0;

err = sum(bsxfun(@and,~bsxfun(@eq,T,P),cvblocks)) ./ sum(cvblocks);

hit = sum(bsxfun(@and,bsxfun(@and,T,P),cvblocks));
fa = sum(bsxfun(@and,bsxfun(@and,~T,P),cvblocks));

hitr = hit./sum(bsxfun(@and,T,cvblocks));
far = fa./sum(bsxfun(@and,~T,cvblocks));
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

counts.targets = sum(bsxfun(@and,T,cvblocks));
counts.distractors = sum(bsxfun(@and,~T,cvblocks));
counts.hits = hit;
counts.hitrate = hitr;
counts.falsealarms = fa;
counts.falsealarmrate = far;
