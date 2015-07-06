function M = group2mat(G)
    n = cellfun('prodofsize', G(:));
    s = [0;cumsum(n)];
    M = nan(numel(G), max(n));
    for i = 1:numel(G)
        a = n(i);
        b = s(i) + 1;
        c = n(i) + s(i);
        M(i,1:a) = b:c;
    end
end