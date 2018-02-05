function y = classifier_error(x,p)
% CLASSIFIER_ERROR Error function for logistic classifiers
%
% INPUT
% x : the target matrix
% p : the predicted matrix
%
% OUTPUT
% y : error term (the loss value)
%
% NOTE
% This is not the logistic loss function. What this returns is
% misclassification error, adjusted for unequal sample sizes.
    pp = p > 0;
    pn = p <= 0;
    xp = x > 0;
    xn = x <= 0;

    poscount = nnz(xp);
    negcount = nnz(xn);

    truepos = xp & pp;
    falsepos = xn & pp;
    trueneg = xn & pn;
    falseneg = xp & pn;

    hitrate = nnz(truepos) / poscount;
    rejectrate = nnz(trueneg) / negcount;

    accuracy = (hitrate + rejectrate) / 2;

    y = 1 - accuracy;
end
