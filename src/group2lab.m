function labs = group2lab(G)
    n = cellfun('prodofsize', G(:));
    M = zeros(max(n),numel(G));
    for i = 1:numel(G)
        a = n(i);
        M(1:a,i) = ones(a,1);
    end
    M = bsxfun(@times, M, 1:numel(G));
    labs = M(M>0);
end
