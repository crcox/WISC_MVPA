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
    %   Author: Chris Cox, 7/8/2015
    % Older releases of matlab do not have the repelem() builtin. If that is the case, use my crappy version instead.
    mlinfo = ver('MATLAB');
    mlversion = sscanf(mlinfo.Version,'%d.%d');
    if (mlversion(1) <= 8) && (mlversion(2) < 5)
        repelem_ = @repelem_crc;
    else
        repelem_ = @repelem;
    end

    p = inputParser();
    addRequired(p, 'coords',   @isnumericlike);
    addRequired(p, 'diameter', @isnumeric);
    addRequired(p, 'overlap',  @isnumeric);
    addRequired(p, 'shape',    @isvalidshape);
    parse(p, coords, diameter, overlap, shape);

    coords   = p.Results.coords;
    diameter = p.Results.diameter(:)';
    overlap  = p.Results.overlap(:)';
    shape    = p.Results.shape;
    if iscell(coords)
        d = cellfun(@(x) size(x,2), coords);
        n = cellfun(@(x) size(x,1), coords);
        n = n(:)'; % force row vec
        N = length(coords(:));
        assert(all(d(1) == d));
        coords = cell2mat(coords(:));
        sid = repelem_(1:N, n);
        ncum = uint32([0, cumsum(n)]);
    else
        [n,d] = size(coords);
        N = 1;
    end

    assert(length(diameter) == 1 || length(diameter) == d);
    assert(length(overlap) == 1 || length(overlap) == d);
    assert(all(overlap <= diameter));

    r = diameter ./ 2;
    x = diameter - overlap;

    radii = true(size(coords,1),1);
    G = cell(size(coords,1),N);
    if any(r == 0)
        if N > 1
            for j = 1:N
                a = ncum(j) + 1;
                b = ncum(j+1);
                G(a:b, j) = num2cell(uint32(1:n(j)));
                G(a:b, j~=[1:N]) = {uint32(zeros(0,1))};
            end
        else
            G(:, 1) = num2cell(uint32(1:n(1)));
        end

    else
        for i = 1:size(coords,1);
            if radii(i) == true;
                switch lower(shape)
                    case 'sphere'
                        tmp = bsxfun(@minus,coords, coords(i,:)).^2;
                        z = sum(bsxfun(@rdivide, tmp, r.^2),2) < 1;
                        g = uint32(find(z));
                        if N > 1
                            s = sid(z);
                            for j = 1:N
                                G{i,j} = g(s==j) - ncum(j);
                            end
                        else
                            G{i,1} = g;
                        end
                        z = sum(bsxfun(@rdivide, tmp, x.^2),2) < 1;
                        radii(z) = false;
                    case 'cube'
                        cmax = coords(i,:) + r;
                        cmin = coords(i,:) - r;
                        z = all(bsxfun(@lt,coords,cmax) & bsxfun(@gt,coords,cmin),2);
                        g = uint32(find(z));
                        if N > 1
                            s = sid(z);
                            for j = 1:N
                                G{i,j} = g(s==j) - ncum(j);
                            end
                        else
                            G{i,1} = g;
                        end
                        cmax = coords(i,:) + x;
                        cmin = coords(i,:) - x;
                        z = all(bsxfun(@lt,coords,cmax) & bsxfun(@gt,coords,cmin),2);
                        radii(z) = false;
                    case 'unitary'
                        warning('coordGrouping:ignoringArgs','The unitary shape implies diameter = overlap = 0. Nonzero values provided are being ignored.');
                        if N > 1
                            for j = 1:N
                                a = ncum(j) + 1;
                                b = ncum(j+1);
                                G(a:b, j) = num2cell(uint32(1:n(j)));
                                G(a:b, j~=[1:N]) = {uint32(zeros(0,1))};
                            end
                        else
                            G(:, 1) = num2cell(uint32(1:n(1)));
                        end
                        break

                end
            end
        end
    end
    G(all(cellfun('isempty',G),2),:) = [];
end

function b = isvalidshape(x)
    b = ischar(x) && any(strcmpi(x,{'sphere','cube','unitary'}));
end

function b = isnumericlike(x)
    if iscell(x)
        b = all(cellfun(@isnumeric, x));
    else
        b = isnumeric(x);
    end
end
