function G = coordGrouping(coords, diameter, overlap, shape)
    % COORDGROUPING  Form groups of coordinates based on definitions.
    %   G = coordGrouping(coords, diameter, overlap, shape) generates groups.
    %   INPUTS
    %   coordinates : a n-by-d matrix of coordinates, or a cell array of such
    %                 matrices. If a cell array is provided, groups are formed
    %                 with respect to the complete set of coordinates, across
    %                 all cells, and then sorted out by coordinate cell. This
    %                 ensures that
    %   diameter    : the distance across shape along each axis. If only one
    %                 diameter is provided, it will apply to all d axes. If
    %                 diameter is zero, groups will each include a single
    %                 voxel.
    %   overlap     : the amount that each group should overlap.
    %   shape       : can be either "cube" or "sphere", although technically
    %                 elipsoids and cuboids are possible by specifying
    %                 different diameters along each dimension.
    %   OUTPUTS
    %   G : A matrix of cells, where rows correspond to groups and there is a
    %       column for each cell of coordinates provided. Each cell contains
    %       indexes that point to rows in the corresponding cell of coords.
    %
    %   See also ndmovingwindow, ndcoord.
    %   Author: Chris Cox, 20 April 2018

    p = inputParser();
    addRequired(p, 'coords',   @(x) iscell(x) && isnumeric(x{1}));
    addRequired(p, 'diameter', @isnumeric);
    addRequired(p, 'overlap',  @isnumeric);
    addRequired(p, 'shape',    @(x) any(strcmpi(x,{'cube','sphere','unitary'})));
    parse(p, coords, diameter, overlap, shape);

    coords   = p.Results.coords;
    diameter = p.Results.diameter(:)';
    overlap  = p.Results.overlap(:)';
    shape    = p.Results.shape;
    
    d = cellfun(@(x) size(x,2), coords);
    n = cellfun(@(x) size(x,1), coords);
    n = n(:)'; % force row vec
    N = length(coords(:));
    assert(all(d(1) == d));
    ncum = uint32([0, cumsum(n)]);
    
    assert(length(diameter) == 1 || all(length(diameter) == d));
    assert(length(overlap) == 1 || all(length(overlap) == d));
    assert(all(overlap <= diameter));

    d = d(1);
    if numel(diameter) == 1
        diameter = repmat(diameter, d, 1);
    end
    
    x = (diameter - overlap);% ./ 2;
    r = diameter ./ 2;

    extents.min = floor(min(cell2mat(coords)));
    extents.max = ceil(max(cell2mat(coords)));
    g = cell(d,1);
    for i = 1:d
        a = extents.min(i);
        b = extents.max(i);
        % g{i} will contain the onset of each group in a particular
        % dimension.
        g{i} = a:x(i):b;
        if strcmpi(shape, 'sphere')
            % If we are generating spheres, the points we are ultimately
            % generating are going to be the radii. If this point is
            % further than r(i) away from b, then we should add an
            % additional point to compose an additional radii.
            if (g{i}(end) + r(i)) < b
                g{i} = [g{i},g{i}(end)+x(i)];
            end
        end
    end
    gg = cell(d,1);
    [gg{:}] = ndgrid(g{:});
    gg = cellfun(@(x) x(:), gg, 'UniformOutput', 0);
    gg = cat(2,gg{:});
    
    G = cell(size(gg,1), N);
    switch lower(shape)
        case 'sphere'

            for j = 1:N
                for i = 1:size(gg,1)
                    di = distance_to_reference(gg(i,:), coords{j});
                    z = all(bsxfun(@le,di(:), r(:)'), 2);
                    G{i,j} = uint32(find(z));        
                end
            end

        case 'cube'
            for j = 1:N
                gd = bsxfun(@plus, gg, diameter(:)');
                for i = 1:size(gg,1)
                    z = bsxfun(@ge, coords{j}, gg(i,:)) & bsxfun(@lt, coords{j}, gd(i,:));
                    z = all(z, 2);
                    G{i,j} = uint32(find(z));        
                end
            end
            
        case 'unitary'
            warning('coordGrouping:ignoringArgs','The unitary shape implies diameter = overlap = 0. Nonzero values provided are being ignored.');
            if N > 1
                for j = 1:N
                    a = ncum(j) + 1;
                    b = ncum(j+1);
                    G(a:b, j) = num2cell(uint32(1:n(j)));
                    G(a:b, j~=(1:N)) = {uint32(zeros(0,1))};
                end
            else
                G(:, 1) = num2cell(uint32(1:n(1)));
            end
    end
    G(all(cellfun('isempty',G),2),:) = [];
end

function d = distance_to_reference(reference, points)
    z = zeros(size(points,1),1);
    for i = 1:numel(reference)
        xp=reference(i);
        x = points(:,i);
        y = ((xp-x).^2);
        z = z + y;
    end
    d = sqrt(z);
end