function [s,ix] = selectbyfield(S, varargin)
% SELECTBYFIELD Subset structured array matching on fields within array
% USAGE
% Example Data:
%   S = struct('label', {'A','A','B','B'}, 'dim', {1,2,1,2});
% Use cases:
%   s = selectbyfield(S, 'label', 'A');
%   s = selectbyfield(S, 'label', {'B'});
%   s = selectbyfield(S, 'label', 'A', 'dim', 1);
%   s = selectbyfield(S, 'label', {'A','B'}, 'dim', 2);
%   s = selectbyfield(S, 'label', 'A', 'dim', [1,2]);

    args = reshape(varargin, 2, []);
    z = true(1, numel(S));
    for i = 1:size(args,2)
        field = args{1, i};
        value = args{2, i};
        if iscell(value) || (~ischar(value) && numel(value) > 1)
            if ~iscell(value)
                value = num2cell(value);
            end
            y = false(1, numel(S));
            for j = 1:numel(value)
                y = y | lookup_value_in_cellarray(value{j}, {S.(field)});
            end
        else
            y = lookup_value_in_cellarray(value, {S.(field)});
        end
        z = z & y;
    end
    s = S(z);
    ix = find(z);
end

function y = lookup_value_in_cellarray(value, cellarray)
    if isnumeric(value)
        y = cellfun(@(x) isequal(value, x), cellarray);
    else
        y = strcmp(value, cellarray);
    end
end
