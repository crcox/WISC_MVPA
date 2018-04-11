function SL = generate_searchlights_0(xyz,radius)
    D = squareform(pdist(xyz));
    SL = cell(size(D,1),1);
    for j = 1:size(D,1)
        SL{j} = find(D(j,:) <= radius);
    end
end
