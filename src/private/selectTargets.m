function Y = selectTargets(metadata, type, label, source, metric, rowfilter)
    Y = cell(1, numel(metadata));
    if ~iscell(rowfilter) && numel(metadata) == 1
        rowfilter = {rowfilter};
    end
    for i = 1:numel(metadata);
        t = metadata(i).targets;
        if strcmpi(type, 'similarity')
            z = all([strcmpi(type,{t.type});strcmpi(label,{t.label});strcmpi(source,{t.sim_source});strcmpi(metric,{t.sim_metric})]);
            if any(z)
                Y{i} = metadata(i).targets(z).target(rowfilter{i}, rowfilter{i});
            end
        else
            z = all([strcmpi(type,{t.type});strcmpi(label,{t.label})]);
            if any(z)
                y = metadata(i).targets(z).target(rowfilter{i});
                Y{i} = y(:);
            end
        end
    end
end
