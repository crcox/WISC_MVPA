function [alpha,lambda] = independent2ratio(lamSOS, lamL1)
    alpha = lamSOS / (lamSOS + lamL1);
    lambda = lamSOS + lamL1;
end
