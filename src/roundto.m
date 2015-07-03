function y = roundto(x,multiple)
    y = bsxfun(@rdivide,x,multiple);
    y = bsxfun(@times, round(y), multiple);
    check = unique(y,'rows');
    if length(check) < length(y);
        warning('roundto resulted in %d duplicate records.', length(y) - length(check))
    end
end
