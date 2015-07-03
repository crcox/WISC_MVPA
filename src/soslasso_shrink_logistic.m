%%% shrinkage on product norm

function y = soslasso_shrink_logistic(X,G,groups,lam,alpha)

    GroupSparse = (alpha)*lam;
    OverallSparse = (1-alpha)*lam;
    % step 1: perform soft thresholding
    X_soft = sign(X).*max(abs(X) - OverallSparse,0);

    %step 2: perform group soft thresholding
		if GroupSparse > 0
			M = size(G,1); % number of groups
			X_soft = [X_soft; zeros(1,size(X_soft,2))]; % for the dummy
			G(isnan(G)) = size(X_soft,2);
			Xtemp = sum(X_soft.^2,2); %xtemp is now a vector
			Xtemp = sum(Xtemp(G),2);
			Xtemp = sqrt(Xtemp);
			Xtemp = max(Xtemp - GroupSparse,0); % this is the multiplying factor
			z = Xtemp > 0; % If GroupSparse = 0 and Xtemp(i) = 0, NaN results. 
										 % So, only shrink Xtemp(i) if not already 0.
			Xtemp(z) = Xtemp(z)./(Xtemp(z) + GroupSparse);
			if (size(Xtemp,1)~=M)
					error('something weird is happening with the group shrinkage \n');
			end
			Xtemp = Xtemp(groups);
			Xtemp = repmat(Xtemp,1,size(X,2));
			y = X_soft(1:end-1,:).*Xtemp;
		else
			y = X_soft;
		end
end
