function [s,ix] = selectbyfield(S, varargin)
% SELECTBYFIELD Subset structured array matching on fields within array
% USAGE
% S = struct('label', {'A','A','B','B'}, 'dim', {1,2,1,2});
% s = selectbyfield(S, 'label', 'A');
% s = selectbyfield(S, 'label', 'A', 'dim', 1);
    args = reshape(varargin, 2, []);
    z = true(1, numel(S));
    for i = 1:size(args,2)
        field = args{1, i};
        value = args{2, i};
        if isnumeric(value)
            z = z & cellfun(@(x) isequal(value, x), {S.(field)});
        else
            z = z & strcmp(value, {S.(field)});
        end
    end
    s = S(z);
    ix = find(z);
end
