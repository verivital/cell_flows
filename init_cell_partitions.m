if opt_randic
    cellSize = 25;
    randsize = 20;
    numpts = 20;    % number of points to triangulate
    pts = cellSize*rand(numpts,m); % broken most of the time
    %pts = ceil(randsize*rand(numpts,m)); % works most of the time
else
    if opt_tesselation_type == OPT_TESSELATION_TRI && opt_fixedic == 9 % figure used in paper
        pts = [1 1; 1     3; 1     3; 1     5; 2     1; 2     5; 3     2; 3     3; 3     4; 4     2; 4     4; 4     5; 5     1; 5     3; 5     4;     5     5];
    else
        %pts = [0 0; 0 1; 1 0; 1 1; 2 0; 0 2; 2 2; 1 2; 2 1];
        %pts = [0 0; 0 1; 1 0; 1 1; 2 0; 0 2; 1 2; 2 1; 2 2; 1 3; 3 1; 3 3; 0 3; 3 0; 3 2; 2 3];
        %pts = [0 0; 0 1; 1 0; 1 1; 2 0; 0 2; 1 2; 2 1; 2 2; 1 3; 3 1; 3 3; 0 3; 3 0; 3 2; 2 3; 1 4; 4 1; 3 4; 4 3; 4 4; 0 4; 4 0; 4 2; 2 4];
        % % generate evenly spaced grid
        % row = [0:6]';
        % pts = [];
        % for i = 1 : length(row)
        %     col = circshift(row,i);
        %     pts = [pts; row, col];
        % end

        % generating regular grid
        lowest = 0;
        highest = 5;
        step = 1;
        [x, y] = meshgrid(lowest:step:highest,lowest:step:highest);
        k = 1;
        num = ceil( (highest - lowest) / step);
        for i = 1 : num + 1
            for j = 1 : num + 1
                pts(k,1) = x(i,j);
                pts(k,2) = y(i,j);
                k = k + 1;
            end
        end
    end
end

% ensure points are unique
pts = unique(pts, 'rows');
pts = sortrows(pts);

if opt_triangulation
    % todo: add constraints describing polygon being triangularized
    dt = DelaunayTri(pts(:,1),pts(:,2));
    dtTriangulationSorted = sortrows(dt.Triangulation);
    inside = dt.inOutStatus();
    % Construct a TriRep object to represent the domain triangles
    tr = TriRep(dt(inside, :), sortrows(dt.X));

    % Construct a set of edges that join the circumcenters of neighboring triangles; the additional logic constructs a unique set of such edges.
    N = size(tr,1); % number of cells/triangles
    diam = N/2;
else
    switch opt_tesselation_type
        case OPT_TESSELATION_TRI
            % generate uniform grid
            N = ((2*num)^2)/2; % num is the number of squares between low and high; N is number of total cells
            diam = 2*ceil( sqrt(3)/2 * sqrt(N)) + length(failed); % rough estimate
            numadd = num + 1;
            j = 1;
            for i = 1 : N
                if mod(i,2) ~= 0 % odd 
                    Cell(i).vertices = [pts(j,:); pts(j + 1,:); pts(j + numadd,:)];
                else % even
                    Cell(i).vertices = [pts(j + numadd,:); pts(j + 1,:); pts(j + 1 + numadd,:)];
                    j = j + 1;
                end

                % increment j twice when we wrap around
                if mod(i,2*(num)) == 0
                    j = j + 1;
                end

                % break out if we will exceed number of points we have
                % note that we iterate over N, but there are generally fewer points
                % than cells (we use the same points for several cells)
                if j + numadd + 1 > length(pts)
                    break;
                end
            end
        case OPT_TESSELATION_EQTRI_LINE
            primitive_cells = 1;
            copies = max(targets);
            N = copies * primitive_cells;
            diam = 2*N; % line graph, worst case

            baseX = 0;
            baseY = 0;
            side = 1;
            Cell(1).vertices = [baseX, baseY; baseX + side, baseY; baseX + side/2, baseY + cosd(60)*side];
            Cell(1).num_sides = 3;
            Cell(1).ptype = OPT_TESSELATION_TRI;

            c = 1;
            numDir = 1;
            cx = 0;
            cv0 = [side/2, 0]; % x only
            for i = primitive_cells + 1 : N
                ci = mod(i-1,primitive_cells)+1;
                Cell(i).vertices = Cell(ci).vertices + ones(Cell(ci).num_sides,1)*(i-1)*cv0;
                
                % flip
                if mod(i,2) == 0
                    Cell(i).vertices(1,2) = 0.5;
                    Cell(i).vertices(2,2) = 0.5;
                    Cell(i).vertices(3,2) = 0;
                end
                Cell(i).num_sides = Cell(ci).num_sides;
                Cell(i).ptype = Cell(ci).ptype;
            end
        % TODO:
        case OPT_TESSELATION_RTRI_LINE
            primitive_cells = 1;
            copies = max(targets);
            N = copies * primitive_cells;
            diam = 2*N; % line graph, worst case

            baseX = 0;
            baseY = 0;
            side = 1;
            Cell(1).vertices = [baseX, baseY; baseX + side/2, baseY; baseX, baseY + side/2];
            Cell(1).num_sides = 3;
            Cell(1).ptype = OPT_TESSELATION_TRI;

            c = 1;
            numDir = 1;
            cx = 0;
            for i = primitive_cells + 1 : N
                ci = mod(i-1,primitive_cells)+1;
                
                %if mod(i,2) == 0
                %    cv0 = [0, 0]; % x only
                %else
                    cv0 = [side/2, 0]; % x only
                %end
                
                Cell(i).vertices = Cell(ci).vertices;
                
                % flip
                if mod(i,2) == 0
                    theta = 45;
                    Arot = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
                    Cell(i).vertices(1,:) = Arot*Cell(i).vertices(1,:)';
                    Cell(i).vertices(2,:) = Arot*Cell(i).vertices(2,:)';
                    Cell(i).vertices(3,:) = Arot*Cell(i).vertices(3,:)';
                    %Cell(i).vertices(1,2) = 0.5;
                    %Cell(i).vertices(2,2) = 0.5;
                    %Cell(i).vertices(3,2) = 0;
                end
                
                Cell(i).vertices = Cell(i).vertices + ones(Cell(ci).num_sides,1)*(i-1)*cv0;
                
                Cell(i).num_sides = Cell(ci).num_sides;
                Cell(i).ptype = Cell(ci).ptype;
            end
        case OPT_TESSELATION_RECTANGULAR
        case OPT_TESSELATION_SQUARE
            % generate uniform grid
            N = ((num)^2); % num is the number of squares between low and high; N is number of total cells
            diam = 2*ceil( sqrt(N) ) + length(failed); % rough estimate
            numadd = num + 1;
            j = 1;
            for i = 1 : N
                Cell(i).vertices = [pts(j,:); pts(j + 1,:); pts(j + numadd,:); pts(j + numadd + 1,:)];
                j = j + 1;

                % increment j twice when we wrap around
                if mod(i,num) == 0
                    j = j + 1;
                end

                % break out if we will exceed number of points we have
                % note that we iterate over N, but there are generally fewer points
                % than cells (we use the same points for several cells)
                if j + numadd + 1 > length(pts)
                    break;
                end
            end
        case OPT_TESSELATION_SQUARE_LINE
            primitive_cells = 1;
            copies = max(targets);
            N = copies * primitive_cells;
            diam = 2*N; % line graph, worst case

            baseX = 0;
            baseY = 0;
            sideX = 1;
            sideY = 1;
            Cell(1).vertices = [baseX, baseY; baseX, baseY + sideY; baseX + sideX, baseY; baseX + sideX, baseY + sideY];
            Cell(1).num_sides = 4;
            Cell(1).ptype = OPT_TESSELATION_SQUARE;

            c = 1;
            numDir = 1;
            cx = 0;
            cv0 = [sideX, 0]; % x only
            for i = primitive_cells + 1 : N
                ci = mod(i-1,primitive_cells)+1;
                Cell(i).vertices = Cell(ci).vertices + ones(Cell(ci).num_sides,1)*(i-1)*cv0;
                Cell(i).num_sides = Cell(ci).num_sides;
                Cell(i).ptype = Cell(ci).ptype;
            end
        case OPT_TESSELATION_SNUB_SQUARE_TILING
            primitive_cells = 6; % primitive tile has 6 cells
            copies = 4;
            N = copies * primitive_cells;
            diam = 2*ceil( sqrt(3)/2* sqrt(N) ) + length(failed); % rough estimate

            basePoint = 100;
            side = 2;
            angledPoint = side*sqrt(3)/2;
            translationPoint = side*1/2*(1 + sqrt(3));
            Cell(1).vertices = [basePoint, 0; basePoint-(side/2), angledPoint; basePoint+(side/2), angledPoint];
            Cell(1).num_sides = 3;
            Cell(1).ptype = OPT_TESSELATION_TRI;
            Cell(2).vertices = [basePoint-(side/2), angledPoint; basePoint+(side/2), angledPoint; basePoint, angledPoint*2];
            Cell(2).num_sides = 3;
            Cell(2).ptype = OPT_TESSELATION_TRI;
            Cell(3).vertices = [100-(side/2), angledPoint; basePoint - translationPoint, translationPoint; 100, angledPoint*2; basePoint - angledPoint, 2*angledPoint + side/2];
            Cell(3).num_sides = 4;
            Cell(3).ptype = OPT_TESSELATION_SQUARE;
            Cell(4).vertices = [100+(side/2), angledPoint; basePoint + translationPoint, translationPoint; 100, angledPoint*2; basePoint + angledPoint, 2*angledPoint + side/2];
            Cell(4).num_sides = 4;
            Cell(4).ptype = OPT_TESSELATION_SQUARE;
            Cell(5).vertices = [basePoint, angledPoint*2; basePoint, angledPoint*2 + side; basePoint - angledPoint, 2*angledPoint + side/2];
            Cell(5).num_sides = 3;
            Cell(5).ptype = OPT_TESSELATION_TRI;
            Cell(6).vertices = [basePoint, angledPoint*2; basePoint, angledPoint*2 + side; basePoint + angledPoint, 2*angledPoint + side/2];
            Cell(6).num_sides = 3;
            Cell(6).ptype = OPT_TESSELATION_TRI;

            c = 1;
            cv = [translationPoint,translationPoint];
            for i = primitive_cells + 1 : N
                if mod(i-1,2*primitive_cells) == 0
                    cv(1) = -cv(1); % add and subtract
                end
                % todo: next part is fairly specific for 4 copies
                if mod(i-1,3*primitive_cells) == 0
                    c = c + 1;
                    cv(1) = 0;
                end

                ci = mod(i-1,primitive_cells)+1;
                Cell(i).vertices = Cell(ci).vertices + ones(Cell(ci).num_sides,1)*c*cv;
                Cell(i).num_sides = Cell(ci).num_sides;
                Cell(i).ptype = Cell(ci).ptype;
            end
        case OPT_TESSELATION_PARALLELOGRAM
            primitive_cells = 1; % primitive tile has 6 cells
            copies = 16;
            N = copies * primitive_cells;
            diam = 2*ceil( sqrt(N) ) + length(failed); % rough estimate

            baseX = 0;
            baseY = 0;
            sideX = 2;
            sideY = 2;
            angle = 30; %todo: used as a global, make clearer
            angleX = sind(angle)*sideX;
            angleY = cosd(angle)*sideY;
            Cell(1).vertices = [baseX, baseY; baseX + angleX, baseY + angleY; baseX + sideX, baseY; baseX + sideX + angleX, baseY + angleY];
            Cell(1).num_sides = 4;
            Cell(1).ptype = OPT_TESSELATION_PARALLELOGRAM;

            c = 1;
            numDir = 2;
            cx = 0;
            cy = 0;
            cv0 = [sideX, 0]; % x only
            cv1 = [angleX, angleY]; % y only
            square = floor(sqrt(N));
            for i = primitive_cells + 1 : N
                cv = cv0*mod(i-1,square) + cy*cv1;
                
                if mod(i,square) == 0
                    cy = cy + 1;
                end

                ci = mod(i-1,primitive_cells)+1;
                Cell(i).vertices = Cell(ci).vertices + ones(Cell(ci).num_sides,1)*cv;
                Cell(i).num_sides = Cell(ci).num_sides;
                Cell(i).ptype = Cell(ci).ptype;
            end
        case OPT_TESSELATION_HEXAGON
            primitive_cells = 1; % primitive tile has 6 cells
            copies = 16;
            N = copies * primitive_cells;
            diam = 2*ceil( sqrt(N) ) + length(failed); % rough estimate

            side = 2;
            angle = 360/6; %todo: used as a global, make clearer
            angleX = sind(angle)*side;
            angleY = cosd(angle)*side;
            baseX = 1/2*side + angleY;
            baseY = sind(angle) * side;
            %baseX = 0;
            %baseY = 0;
            
            Cell(1).vertices = zeros(6,2);
            for sd = 1 : 6
                Cell(1).vertices(sd,:) = [baseX + side * cosd(sd * angle), baseY + side * sind(sd * angle)];
            end
            Cell(1).vertices
            Cell(1).num_sides = 6;
            Cell(1).ptype = OPT_TESSELATION_HEXAGON;

            c = 1;
            numDir = 2;
            cx = 0;
            cy = 0;
            cv0 = [3/2*side, baseY]; % x only
            cv1 = [0, 2*baseY]; % y only
            %cv2 = [
            square = floor(sqrt(N));
            for i = primitive_cells + 1 : N
                cv = cv0*mod(i-1,square) + cy*cv1;
                
                if mod(i,square) == 0
                    cy = cy + 1;
                end

                ci = mod(i-1,primitive_cells)+1;
                Cell(i).vertices = Cell(ci).vertices + ones(Cell(ci).num_sides,1)*cv;
                Cell(i).num_sides = Cell(ci).num_sides;
                Cell(i).ptype = Cell(ci).ptype;
            end
        case OPT_TESSELATION_HEXAGON_DIAMOND % todo: add 2nd cell type
            primitive_cells = 1; % primitive tile has 6 cells
            copies = 12;
            N = copies * primitive_cells;

            side = 2;
            angle = 360/6; %todo: used as a global, make clearer
            angleX = sind(angle)*side;
            angleY = cosd(angle)*side;
            baseX = 1/2*side + angleY;
            baseY = sind(angle) * side;
            %baseX = 0;
            %baseY = 0;
            
            Cell(1).vertices = zeros(6,2);
            for sd = 1 : 6
                Cell(1).vertices(sd,:) = [baseX + side * cosd(sd * angle), baseY + side * sind(sd * angle)];
            end
            Cell(1).vertices
            Cell(1).num_sides = 6;
            Cell(1).ptype = OPT_TESSELATION_HEXAGON;

            c = 1;
            numDir = 2;
            cx = 0;
            cy = 0;
            cv0 = [2*baseX, 0]; % x only
            cv1 = [0, 2*baseY]; % y only
            square = floor(sqrt(N));
            for i = primitive_cells + 1 : N
                cv = cv0*mod(i-1,square) + cy*cv1;
                
                if mod(i,square) == 0
                    cy = cy + 1;
                end

                ci = mod(i-1,primitive_cells)+1;
                Cell(i).vertices = Cell(ci).vertices + ones(Cell(ci).num_sides,1)*cv;
                Cell(i).num_sides = Cell(ci).num_sides;
                Cell(i).ptype = Cell(ci).ptype;
            end
    end
end