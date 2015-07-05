function val = nonsmooth_eval(Wz,lam,alpha,groups)
	GroupSparse = (alpha)*lam;
	OverallSparse = (1-alpha)*lam;

	% L1 Norm (Lasso component)
	val = sum(sum(abs(Wz)))*OverallSparse;

	% Group sparse component
  if max(max(groups)) > size(Wz,1)
    Wz = [Wz; zeros(1,size(Wz,2))]; % for the dummy
  end
	groups(isnan(groups)) = size(Wz,1);
	Wz = sum(Wz.^2,2); 
	Wz = sum(Wz(groups),2);
	Wz = sqrt(Wz);
	Wz = sum(Wz)*GroupSparse;

	% Combine them
	val = val + Wz;
end