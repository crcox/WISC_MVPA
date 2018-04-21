function [ ModelInstances ] = SOSGroupLoad( ModelInstances )
%SOSGROUPLOAD Summary of this function goes here
%   Detailed explanation goes here
    [MIT,ia,ib] = unique(table( ...
        {ModelInstances.sosgroups}', ...
        'VariableNames', {'sosgroups'}));
    iA = ia(ib);
    G = cell(size(MIT,1),1);
    for ii = 1:size(MIT,1)
        tmp = load(MIT.sosgroups{ii},'G');
        G{ii} = tmp.G;
    end
    for ii = 1:numel(G)
        z = iA == ii;
        [ModelInstances(z).G] = deal(G{ii});
    end

end
