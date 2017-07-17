function v = group2vec(G)
    % NB: This works because (:) notation will always produce a column
    % vector, so cell2mat is sure to produce a column vector, as well.
    G = cellfun(@(x) x(:), G(:), 'Unif', 0);
    v = cell2mat(G);
end
