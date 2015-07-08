xyz = ndcoord(1:30,1:30,1:60);
xyz = bsxfun(@minus, xyz, mean(xyz));
d = pdist2(xyz,[0,0,0]);
xyz = xyz(d<15,:);

groupdiameter = 6;
r = groupdiameter / 2;
overlap = 3;
x = groupdiameter - overlap;
radii = true(size(xyz,1),1);
G = cell(size(xyz,1),1);
shape = 'cube';
for i = 1:size(xyz,1);
    if radii(i) == true;
        switch shape
            case 'sphere'
                d = pdist2(xyz,xyz(i,:));
                G{i} = find(d < r);
                radii(d < x) = false;
            case 'cube'
                xyzmax = xyz(i,:) + r;
                xyzmin = xyz(i,:) - r;
                z = all(bsxfun(@lt,xyz,xyzmax) & bsxfun(@gt,xyz,xyzmin),2);
                G{i} = find(z);
                xyzmax = xyz(i,:) + x;
                xyzmin = xyz(i,:) - x;
                z = all(bsxfun(@lt,xyz,xyzmax) & bsxfun(@gt,xyz,xyzmin),2);
                radii(z) = false;
        end
    end
end
G(cellfun('isempty',G)) = [];
    
D = pdist(xyz);

D