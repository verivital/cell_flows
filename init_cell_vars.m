% Iterate over all N cells, initializing an array of structs of length N
%   * generate the neighbor graph for the cells using data structures generated
%     by the DelaunayTri command
%   * Note: cell is a keyword in Matlab, so make sure you use the word Cell to
%     when using this variable
for i = 1 : N
    if opt_triangulation
        %Cell(i).nbrs = tr.neighbors(i);                 % copy neighbors from triangulation data (seems buggy, use our way below)
        %Cell(i).incenter = ic(i,:);                     % incenter of the cell (for plotting, etc); this is a vector
        %Cell(i).R = icr(i,:);
    end
    Cell(i).nbrs = []; % init
    if opt_tesselation_type == OPT_TESSELATION_TRI || opt_tesselation_type == OPT_TESSELATION_EQTRI_LINE || opt_tesselation_type == OPT_TESSELATION_RTRI_LINE
        Cell(i).num_sides = 3;                      % all triangles
        Cell(i).ptype = OPT_TESSELATION_TRI;
    elseif opt_tesselation_type == OPT_TESSELATION_SQUARE || opt_tesselation_type == OPT_TESSELATION_SQUARE_LINE
        Cell(i).num_sides = 4;                      % all squares
        Cell(i).ptype = OPT_TESSELATION_SQUARE;
    end
    Cell(i).detectNextFailed = 0;                      % round at which a cell's next pointer was detected as being failed---it must wait until: k > detectNextFailed + firstRealRound before resuming operation
    Cell(i).dist = inf * ones(length(targets),1);   % infinite distance to every target initially
    Cell(i).next = i * ones(length(targets),1);     % point at oneself initially
    Cell(i).signal = i;                             % signal used by Cell i to indicate to a neighbor that it is safe to move in the direction of Cell i (i.e., no entities on Cell i along their common border)
    Cell(i).lastSideSignal = 0;                     % last side allowed to move
    Cell(i).token = [];                             % token used to fairly signal each neighbor
    Cell(i).failed = 0;                             % initially cells are not faulty
    Cell(i).Entities = [];                          % set of entities on Cell i
    Cell(i).color = 0; % null color, for removing entities on targets
    Cell(i).etype = [];                             % types of entities on this cell
    Cell(i).lastEtype = [];                         % last entity types on this cell
    Cell(i).throughput = zeros(NT,1);
    for tt = 1 : length(targets)
        %Cell(i).epsPlot(tt,:) = [epsPlot*rand(1,1), epsPlot*rand(1,1)];              % small spacing off the center for easier visualizing (see plotSystem)
        if mod(tt,2) == 0
            Cell(i).epsPlot(tt,:) = [-epsPlot, -epsPlot];              % small spacing off the center for easier visualizing (see plotSystem)
        else
            Cell(i).epsPlot(tt,:) = [epsPlot, epsPlot];              % small spacing off the center for easier visualizing (see plotSystem)
        end
        Cell(i).path(tt).nextTried = [];
        Cell(i).path(tt).lockset = [];
        Cell(i).path(tt).lockcolors = [];
        Cell(i).path(tt).ids = [];
        Cell(i).path(tt).int = [];
        Cell(i).path(tt).nlock = [];    % set of colors which need to lock this cell to move
        Cell(i).path(tt).lock = 0;      % lock for permission to move along a path
        Cell(i).path(tt).rlock = 0;     % bit to reset lock after entity transfer
        Cell(i).path(tt).dlock = 1;     % bit to specify done with this lock
        Cell(i).path(tt).prev = []; % set of things pointing at i
    end
    
    Cell(i).ud = [];
    Cell(i).speed = L * 0.99; % needs to be < L for safety
    
    if opt_triangulation
        Cell(i).vindex = dtTriangulationSorted(i ,:); %the member variable is a 1*3 vector that contains what points the triangle consists of
        Cell(i).vertices = [pts(Cell(i).vindex(1,1), :) ; pts(Cell(i).vindex(1,2), :) ; pts(Cell(i).vindex(1,3), :) ]; % initilize the coordinates with vertices with actual coordinates
    
        Cell(i).side(1).vertices = [pts(Cell(i).vindex(1,1),:); pts(Cell(i).vindex(1,2),:)];
        Cell(i).side(2).vertices = [pts(Cell(i).vindex(1,1),:); pts(Cell(i).vindex(1,3),:)];
        Cell(i).side(3).vertices = [pts(Cell(i).vindex(1,2),:); pts(Cell(i).vindex(1,3),:)];
    else
        if Cell(i).num_sides == 3
            Cell(i).side(1).vertices = [Cell(i).vertices(1,:); Cell(i).vertices(2,:)]; % side 1
            Cell(i).side(2).vertices = [Cell(i).vertices(1,:); Cell(i).vertices(3,:)]; % side 2
            Cell(i).side(3).vertices = [Cell(i).vertices(2,:); Cell(i).vertices(3,:)]; % side 3
        elseif Cell(i).num_sides == 4
            Cell(i).side(1).vertices = [Cell(i).vertices(1,:); Cell(i).vertices(2,:)]; % side 1
            Cell(i).side(2).vertices = [Cell(i).vertices(2,:); Cell(i).vertices(4,:)]; % side 2
            Cell(i).side(3).vertices = [Cell(i).vertices(1,:); Cell(i).vertices(3,:)]; % side 3
            Cell(i).side(4).vertices = [Cell(i).vertices(3,:); Cell(i).vertices(4,:)]; % side 4
        elseif Cell(i).num_sides == 6
            Cell(i).side(1).vertices = [Cell(i).vertices(1,:); Cell(i).vertices(2,:)]; % side 1
            Cell(i).side(2).vertices = [Cell(i).vertices(2,:); Cell(i).vertices(4,:)]; % side 2
            Cell(i).side(3).vertices = [Cell(i).vertices(2,:); Cell(i).vertices(3,:)]; % side 3
            Cell(i).side(4).vertices = [Cell(i).vertices(4,:); Cell(i).vertices(5,:)]; % side 4
            Cell(i).side(5).vertices = [Cell(i).vertices(5,:); Cell(i).vertices(6,:)]; % side 5
            Cell(i).side(6).vertices = [Cell(i).vertices(6,:); Cell(i).vertices(1,:)]; % side 6
        end
    end
    
    throughput = zeros(NT,1);
    
    tmpsum = 0;
    for s = 1 : Cell(i).num_sides
        tmpsum = tmpsum + Cell(i).vertices(s,:);
    end
    Cell(i).centroid = tmpsum ./ Cell(i).num_sides;     % point in the plane of the center of the cell (for plotting, etc); this is a vector

    for s = 1 : Cell(i).num_sides
        Cell(i).side(s).signal = i;
        Cell(i).side(s).boundary = 1; % will reset later if the side is not actually on the boundary
        Cell(i).side(s).midpoint = (Cell(i).side(s).vertices(1,:) + Cell(i).side(s).vertices(2,:)) / 2;
        % computer normal vector pointing out of the cell from side s
        Cell(i).side(s).vectorNormal = [-(Cell(i).side(s).vertices(1,2) - Cell(i).side(s).vertices(2,2)), Cell(i).side(s).vertices(1,1) - Cell(i).side(s).vertices(2,1)];
        Cell(i).side(s).vectorNormal = Cell(i).side(s).vectorNormal / norm(Cell(i).side(s).vectorNormal, 2); % normalize
        
        % make movement vectors in the correct directions (TODO: generalize, this is rather a big hack based on which side corresponds to which vector, etc.)
        if Cell(i).num_sides == 3
            Cell(i).side(s).oppositeVertex = setdiff(Cell(i).vertices, Cell(i).side(s).vertices, 'rows');
            % if necessary, reverse normal vector to ensure it points away from Cell i (it should be farther from the other vertice of triangle T_i)
            if norm( (Cell(i).side(s).midpoint + Cell(i).side(s).vectorNormal) - Cell(i).side(s).oppositeVertex, 2) < norm( Cell(i).side(s).midpoint - Cell(i).side(s).oppositeVertex, 2)
                Cell(i).side(s).vectorNormal = -Cell(i).side(s).vectorNormal;
            end
        elseif Cell(i).num_sides == 4
            if (s == 1 || s == 2) && Cell(i).ptype == OPT_TESSELATION_SQUARE && (opt_tesselation_type == OPT_TESSELATION_SQUARE || opt_tesselation_type == OPT_TESSELATION_SQUARE_LINE)
                Cell(i).side(s).vectorNormal = -Cell(i).side(s).vectorNormal;
            elseif opt_tesselation_type == OPT_TESSELATION_SNUB_SQUARE_TILING
                %Cell(i).ptype == OPT_TESSELATION_SNUB_SQUARE_TILING
                %if (i == 4 || i == 9 || i == 15) && s == 3
                Cell(i).side(s).vectorNormal = -Cell(i).side(s).vectorNormal;
                if mod(i-1,primitive_cells) == 2 && (s == 3 || s == 4)
                    Cell(i).side(s).vectorNormal = -Cell(i).side(s).vectorNormal;
                elseif mod(i-1,primitive_cells) == 3 && (s == 1 || s == 2)
                    Cell(i).side(s).vectorNormal = -Cell(i).side(s).vectorNormal;
                end
            elseif Cell(i).ptype == OPT_TESSELATION_PARALLELOGRAM
                if s == 3
                    Cell(i).side(s).vectorMovement = [-sind(angle), -cosd(angle)];
                elseif s == 2
                    Cell(i).side(s).vectorMovement = [sind(angle), cosd(angle)];
                elseif s == 1
                    Cell(i).side(s).vectorMovement = [-1, 0];
                else
                    Cell(i).side(s).vectorMovement = [1, 0];
                end
            elseif Cell(i).ptype == OPT_TESSELATION_PARALLELOGRAM
                if s == 3
                    Cell(i).side(s).vectorMovement = [-sind(angle), -cosd(angle)];
                elseif s == 2
                    Cell(i).side(s).vectorMovement = [sind(angle), cosd(angle)];
                elseif s == 1
                    Cell(i).side(s).vectorMovement = [-1, 0];
                else
                    Cell(i).side(s).vectorMovement = [1, 0];
                end
            end
            
            if Cell(i).ptype ~= OPT_TESSELATION_PARALLELOGRAM
                Cell(i).side(s).vectorMovement = Cell(i).side(s).vectorNormal;
            end
        else
            Cell(i).side(s).vectorMovement = Cell(i).side(s).vectorNormal;
            
            if s == 6
                Cell(i).side(s).vectorMovement = -Cell(i).side(s).vectorNormal;
            end
        end
        
        Cell(i).side(s).length = norm( Cell(i).side(s).vertices(1,:) - Cell(i).side(s).vertices(2,:), 2); % Euclidean distance
    end
    %~opt_triangulation &&
    % compute coordinates of incenter
    if Cell(i).ptype == OPT_TESSELATION_TRI
        a = Cell(i).side(3).length; % side 3 is opposite of vertex 1 (see line labeled with comment side 1 through side 3 above)
        b = Cell(i).side(2).length; % side 2 is opposite of vertex 2
        c = Cell(i).side(1).length; % side 1 is opposite of vertex 3
        x1 = Cell(i).vertices(1,1);
        x2 = Cell(i).vertices(2,1);
        x3 = Cell(i).vertices(3,1);
        y1 = Cell(i).vertices(1,2);
        y2 = Cell(i).vertices(2,2);
        y3 = Cell(i).vertices(3,2);
        Cell(i).incenter = [(a * x1 + b * x2 + c * x3) / (a + b + c), (a * y1 + b * y2 + c * y3) / (a + b + c)]; % incenter is at the intersection of the three angle bisectors (and is not in general the centroid)
    elseif ~opt_triangulation
        Cell(i).incenter = Cell(i).centroid; % TODO: not in general true, but is true for regular tesselations
    end
    
    if Cell(i).num_sides == 3
        AB = Cell(i).vertices(1,:) - Cell(i).vertices(2,:);
        AC = Cell(i).vertices(1,:) - Cell(i).vertices(3,:);
        BC = Cell(i).vertices(2,:) - Cell(i).vertices(3,:);
        Cell(i).areac = 0.5 * norm(cross([AB,0], [AC,0])); % half of parallelogram method
        Cell(i).area = 0.5 * abs( AB(1) * AC(2) - AB(2) * AC(1)); % since 2d, equals areac
        Cell(i).perimeter = norm(AB,2) + norm(AC,2) + norm(BC,2);
    elseif Cell(i).num_sides == 4
        Cell(i).perimeter = 0;
        longside = -inf;
        shortside = inf;
        if Cell(i).ptype == OPT_TESSELATION_PARALLELOGRAM
            vertOrder = [1 2; 1 3; 3 4; 2 4];
        else
            vertOrder = [1 2; 1 3; 3 4; 2 4];
        end
        for s = 1 : Cell(i).num_sides
            sidelength = norm( Cell(i).vertices(vertOrder(s,1),:) - Cell(i).vertices(vertOrder(s,2),:), 2); % this relies on having an appropriate mapping to find the next adjacent vertex, not just the next vertex
            Cell(i).perimeter = Cell(i).perimeter + sidelength;
            longside = max(longside, sidelength);
            shortside = min(shortside, sidelength);
        end
        
        if Cell(i).ptype == OPT_TESSELATION_SQUARE || Cell(i).ptype == OPT_TESSELATION_RECTANGULAR
            Cell(i).area = shortside * longside; % todo: generalize
        elseif Cell(i).ptype == OPT_TESSELATION_PARALLELOGRAM
            Cell(i).area = longside * shortside * cosd(angle);
        end
    elseif Cell(i).num_sides == 6
        Cell(i).perimeter = 0;
        longside = -inf;
        shortside = inf;
        vertOrder = [1 2; 2 3; 3 4; 4 5; 5 6; 6 1];
        for s = 1 : Cell(i).num_sides
            sidelength = norm( Cell(i).vertices(vertOrder(s,1),:) - Cell(i).vertices(vertOrder(s,2),:), 2); % this relies on having an appropriate mapping to find the next adjacent vertex, not just the next vertex
            Cell(i).perimeter = Cell(i).perimeter + sidelength;
            longside = max(longside, sidelength);
            shortside = min(shortside, sidelength);
        end

        Cell(i).area = 3*sqrt(3)/2*shortside^2; % this assumes regular
    end
    
    [Cell(i).A,Cell(i).b] = vert2con(Cell(i).vertices);
    
    Cell(i).Rp = 2*Cell(i).area / Cell(i).perimeter; % should equal R from triangulation data
    Cell(i).safeRegion.R = Cell(i).Rp - safeDistance;
    
    if opt_tesselation_type == OPT_TESSELATION_HEXAGON
        Cell(i).illegalRegion.R = Cell(i).Rp - 4*L;
    else
        Cell(i).illegalRegion.R = Cell(i).Rp - L;
    end

    Cell(i).illegalRegion.hcoeff = Cell(i).illegalRegion.R / Cell(i).Rp;
    
    for s = 1 : Cell(i).num_sides
        Cell(i).illegalRegion.vertices(s,:) = Cell(i).incenter - ((Cell(i).incenter - Cell(i).vertices(s,:)) * Cell(i).illegalRegion.hcoeff);
        
        % we associate the side as well, used in transfer
        for a = 1 : 2
            Cell(i).illegalRegion.side(s).vertices(a,:) = Cell(i).incenter - ((Cell(i).incenter - Cell(i).side(s).vertices(a,:)) * Cell(i).illegalRegion.hcoeff);
        end
    end
    [Cell(i).illegalRegion.A,Cell(i).illegalRegion.b] = vert2con(Cell(i).illegalRegion.vertices);

    % compute homothety ("shrinking" of the triangle) inside the triangle
    % computation based on this StackExchange post (it had a small bug): http://math.stackexchange.com/questions/17561/how-to-shrink-a-triangle
    % centered from the incenter of the triangle (intersection point of the three angle bisector lines)
    Cell(i).safeRegion.hcoeff = Cell(i).safeRegion.R / Cell(i).Rp;
    if Cell(i).safeRegion.hcoeff > 1 || Cell(i).safeRegion.R <= 0
        'POSSIBLE ERROR: cells are too small when compared to entity radius and safety spacing, change parameters OR only a single entity is allowed per cell, reconfiguring for this case.'
        i
        L
        d
        Cell(i).safeRegion.R
        opt_singleEntity = 1;
        if opt_display
            %call_plotSystem
        end
        %return
    end
    for s = 1 : Cell(i).num_sides
        Cell(i).safeRegion.vertices(s,:) = Cell(i).incenter - ((Cell(i).incenter - Cell(i).vertices(s,:)) * Cell(i).safeRegion.hcoeff);
    end
    
    [Cell(i).safeRegion.A,Cell(i).safeRegion.b] = vert2con(Cell(i).safeRegion.vertices);
end