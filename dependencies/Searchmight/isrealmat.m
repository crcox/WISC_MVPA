% replacement for isrealmat, which some versions of MATLAB do not have

function result = isrealmat(X)

[r,c]= size(X);

result = 1;
result = result & isreal(X);
result = result && isnumeric(X);
result = result && (r>1 | c>1);
