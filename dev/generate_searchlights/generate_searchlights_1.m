function SL = generate_searchlights_1(xyz,radius)
    n = size(xyz,1);
    SL = cell(n,1);
    for j = 1:n
        z = distance_to_reference(xyz(j,:), xyz) <= radius;
        SL{j} = find(z);
    end
end

function d = distance_to_reference(reference, points)
   xp=reference(1); yp=reference(2); zp=reference(3);
   x = points(:,1); y = points(:,2); z = points(:,3);
   d = sqrt(((xp-x).^2) + ((yp-y).^2) + ((zp-z).^2));
end