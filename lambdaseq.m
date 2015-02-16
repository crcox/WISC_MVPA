function seq = lambdaseq(X,Y,alpha,n)
 if iscell(X)
   x = cell2mat(X);
 end

 if iscell(Y)
   y = cell2mat(Y);
 end
%   maxlam = zeros(1,length(X));
%   minlam = zeros(1,length(X));
%   for i = 1:length(X)
%     x = X{i};
%     y = double(Y{i});
%     maxlam(i) = max(x'*y) * alpha;
%     minlam(i) = (1-alpha) * maxlam(i);
%   end
  maxlam = max(x'*y) * alpha;
  minlam = (1-alpha) * maxlam(i);
  seq = log(linspace(exp(min(minlam)),exp(max(maxlam)),10));
end
