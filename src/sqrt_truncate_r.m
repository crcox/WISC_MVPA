function [C, r] = sqrt_truncate_r(S, tau)

%@S: n x n Similarity matrix
%@r: rank tuning parameter
%info: finds square root of S using eigen decompostion and truncates to
%rank r

    %[U,Z] = eig(S);
    [U, Z, ~] = svd(S);

    z = diag(Z);
    n = size(U,2);

    if tau > 1 && (tau - floor(tau)) < 1e-6
        fprintf('Selecting a %d dimensional embedding, treating whole-valued tau > 1 as r.\n',  tau);
        r = tau;
        C = U(:,1:r)*diag(sqrt(z(1:r)));
    else
        fprintf('Picking r such that (norm(S-C*C'',''fro'')/norm(S,''fro'')) <= tau (where tau=%0.2f) ...\n',  tau);
        for r = 1:n
            C = U(:,1:r)*diag(sqrt(z(1:r)));
            objfunc = (norm(S-C*C','fro')/norm(S,'fro'));
            if objfunc <= tau
                break
            end
        end
    end
end
