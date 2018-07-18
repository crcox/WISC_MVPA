function v = repelem_crc(x,n)
    c = [0;cumsum(n(:))];
    v = zeros(c(end),1);
    for i = 1:length(x)
        a = c(i) + 1;
        b = c(i+1);
        v(a:b) = x(i);
    end
end
