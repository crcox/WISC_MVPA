function [lamSOS,lamL1] = ratio2independent(alpha, lambda)
    lamL1 = lambda * (1-alpha);
    lamSOS = lambda * alpha;
end
