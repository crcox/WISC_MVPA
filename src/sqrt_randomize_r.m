function [Crand, r] = sqrt_truncate_r(S, tau)

%@S: n x n Similarity matrix
%@r: rank tuning parameter
%info: finds square root of S using eigen decompostion and truncates to
%rank r

    %[U,Z] = eig(S);
    [U, Z, V] = svd(S);

    z = diag(Z);
    n = size(U,2);

    for r = 1:n
        C = U(:,1:r)*diag(sqrt(z(1:r)));
%         if z(1)>=z(end)
%             C = U(:,1:r)*diag(z(1:r)).^(1/2);
%         else
%             C = U(:,n-r+1:n)*diag(z(n-r+1:n)).^(1/2);
%         end

        %C(abs(C)<0.001)=0;
        objfunc = (norm(S-C*C','fro')/norm(S,'fro'));
  %      disp([r,objfunc]);
        if objfunc <= tau
%            fprintf('Truncated at rank: %d\n',r);
            break
        end
    end
    Crand = randn(size(C)) * Z;
end
