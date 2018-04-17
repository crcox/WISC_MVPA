for i = 1:numel(results)
    p = params(i);
    s = results(i).subject.subject;
    m = metadata(s);
    c = m.coords(2);
    f = selectbyfield(m.filters, 'label', 'colfilter_aud');
    c.ind = c.ind(f.filter);
    c.ijk = c.ijk(f.filter, :);
    c.xyz = c.xyz(f.filter, :);
    results(i).coords = c;
end
    