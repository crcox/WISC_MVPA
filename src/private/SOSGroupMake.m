function [ ModelInstances ] = SOSGroupMake( ModelInstances, xyz )
%SOSGROUPMAKE Summary of this function goes here
%   Detailed explanation goes here
    [diameter_u,~,ib] = unique(cat(1,ModelInstances.diameter),'rows');
    [overlap_u,~,ic] = unique(cat(1,ModelInstances.overlap),'rows');
    [MIT,ia,id] = unique(table( ...
        ib, ...
        ic, ...
        {ModelInstances.shape}', ...
        'VariableNames', {'diameter','overlap','shape'}));
    iA = ia(id);
    MIT.diameter = num2cell(MIT.diameter);
    MIT.overlap = num2cell(MIT.overlap);
    G = cell(size(MIT,1),1);
    for ii = 1:size(MIT,1)
        MIT.diameter{ii} = diameter_u(MIT.diameter{ii},:);
        MIT.overlap{ii} = overlap_u(MIT.overlap{ii},:);
        G{ii} = coordGrouping(xyz, ...
            MIT.diameter{ii}, ...
            MIT.overlap{ii}, ...
            MIT.shape{ii});
    end
    for ii = 1:numel(G)
        z = iA == ii;
        [ModelInstances(z).G] = deal(G{ii});
    end
%     for ii = 1:numel(ModelInstances)
%         diameter = ModelInstances(ii).diameter;
%         overlap = ModelInstances(ii).overlap;
%         shape = ModelInstances(ii).shape;
%         if size(MIT,1) > 1
%             z = all(bsxfun(@eq, cell2mat(MIT.diameter), diameter), 2) & ...
%                 all(bsxfun(@eq, cell2mat(MIT.overlap), overlap), 2) & ...
%                 strcmp(MIT.shape, shape);
%         else
%             z = 1;
%         end 
%         ModelInstances(ii).G = G{z};
%         ModelInstances(ii).subject = [SubjectArray.subject];
%     end

end
