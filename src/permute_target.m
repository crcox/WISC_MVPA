function X = permute_target(Y,method,cvind)
    if ~iscell(Y)
        Y = {Y};
    end
    if nargin < 2
        method = 'simple';
    end 
    if nargin < 3
        cvind = cell(size(Y));
    else
        if ~iscell(cvind)
            cvind = {cvind};
        end
    end
    X = cell(size(Y));
    for i = 1:numel(Y)
        if isempty(cvind{i})
            cvind{i} = ones(size(Y{i}));
        end
        if size(Y{i},1) == 1
            Y{i} = Y{i}';
        end
        cvset = unique(cvind{i})';
        X{i} = zeros(size(Y{i}));
        for j = 1:numel(cvset); % unique(cvind{i})'
            ic = cvset(j);
            y = Y{i}(cvind{i}==ic,:);
            switch method
                case 'simple'
                    X{i}(cvind{i}==ic,:) = simple_perm(y);
                case 'stratified'
                    X{i}(cvind{i}==ic,:) = stratified_perm(y);
            end
        end
    end
end

function x = simple_perm(y)
% Takes y and randomly permutes it
    permix = randperm(size(y,1));
    x = y(permix, :);
end

function x = stratified_perm(y)
% Stratified permuation involves permuting within each category, so that
% each category member is equally likely to be reassigned to each other
% group. Another way to put it is that the permuted result will be as close
% to orthogonal with the true structure as possible, which means that the
% hypothesized structure between the data and this target structure will
% have been eliminated.
    yu = unique(y);
    n = uint16(numel(yu));
    x = zeros(size(y));
    for i = 1:n
        yc = yu(i);
        m = uint16(nnz(y==yc));
        yy = zeros(m,1);
        p = repmat(idivide(m,n),n,1);
        % If the group size is not divisible by the number of groups,
        % select randomly among the group labels to determine the remaining
        % category labels. This way, the inbalance will wash out over
        % multiple permutations.
        r = rem(m,n);
        if r > 0
            ix = randperm(n,r);
            p(ix) = p(ix) + 1;
        end
        cur = 0;
        for j = 1:n
            a = cur + 1;
            b = cur + p(j);
            cur = b;
            yy(a:b) = yu(j);
        end
       x(y==yc) = yy(randperm(m));
    end
end