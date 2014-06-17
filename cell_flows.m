%function target_tp = cell_flows(targets, sources)

OPT_TESSELATION_TRI = 0;
OPT_TESSELATION_SQUARE = 1;
OPT_TESSELATION_SQUARE_LINE = 2;
OPT_TESSELATION_RECTANGULAR = 3;
OPT_TESSELATION_SNUB_SQUARE_TILING = 4;
OPT_TESSELATION_ELONGATED_TRIANGULAR_TILING = 5;
OPT_TESSELATION_PARALLELOGRAM = 6; % parallelogram works, but vectors are different: move entities with a vector in the same direction as the opposite parallel side (e.g., /_/ going up would move with the vector /)
OPT_TESSELATION_HEXAGON = 7;
OPT_TESSELATION_EQTRI_LINE = 8; % line of equilateral triangles
OPT_TESSELATION_RTRI_LINE = 9; % line of right triangles

opt_single = 13; % 0 = batch mode, run all simulations, else value is the number of a simulation below
opt_run = 1; % 0 = batch mode, run all simulations, else value is the number of a run of a simulation below

%MAX_ROUNDS = N + 1;    % needs to be at least N (for routing to have necessarily stabilized), change as necessary
%MAX_ROUNDS = inf;       % run forever

simulations(1).environment = OPT_TESSELATION_SQUARE_LINE;
simulations(1).maxRounds = 2500;
simulations(1).run(1).sources = [1]; % 0 overlap
simulations(1).run(1).targets = [3];

simulations(1).run(2).sources = [1 1]; % 1 overlap
simulations(1).run(2).targets = [3 3];

simulations(1).run(3).sources = [1 1 1]; % 2 overlap
simulations(1).run(3).targets = [3 3 3];

simulations(1).run(4).sources = [1 1 1 1]; % 3 overlap
simulations(1).run(4).targets = [3 3 3 3];

simulations(1).run(5).sources = [1 1 1 1 1]; % 4 overlap
simulations(1).run(5).targets = [3 3 3 3 3];

simulations(1).run(6).sources = [1 1 1 1 1 1]; % 5 overlap
simulations(1).run(6).targets = [3 3 3 3 3 3];

simulations(1).run(7).sources = [1 1 1 1 1 1 1]; % 6 overlap
simulations(1).run(7).targets = [3 3 3 3 3 3 3];

simulations(1).run(8).sources = [1 1 1 1 1 1 1 1]; % 7 overlap
simulations(1).run(8).targets = [3 3 3 3 3 3 3 3];

simulations(1).run(9).sources = [1 1 1 1 1 1 1 1 1]; % 8 overlap
simulations(1).run(9).targets = [3 3 3 3 3 3 3 3 3];


simulations(2).environment = OPT_TESSELATION_SQUARE_LINE;
simulations(2).maxRounds = 2500;
simulations(2).run(1).sources = [8 9]; % square left-right throughput, 0 overlap  (no turns)
simulations(2).run(1).targets = [1 16];

simulations(2).run(2).sources = [8 8]; % square left-right throughput, 1 overlap  (no turns)
simulations(2).run(2).targets = [1 15]; 

simulations(2).run(3).sources = [8 7]; % square left-right throughput, 2 overlap  (no turns)
simulations(2).run(3).targets = [1 14]; 

simulations(2).run(4).sources = [8 6]; % square left-right throughput, 3 overlap  (no turns)
simulations(2).run(4).targets = [1 13];

simulations(2).run(5).sources = [8 5]; % square left-right throughput, 4 overlap  (no turns)
simulations(2).run(5).targets = [1 12];

simulations(2).run(6).sources = [8 4]; % square left-right throughput, 5 overlap  (no turns)
simulations(2).run(6).targets = [1 11];

simulations(2).run(7).sources = [8 3]; % square left-right throughput, 6 overlap  (no turns)
simulations(2).run(7).targets = [1 10];

simulations(2).run(8).sources = [8 2]; % square left-right throughput, 7 overlap  (no turns)
simulations(2).run(8).targets = [1 9];

simulations(2).run(9).sources = [8 1]; % square left-right throughput, 8 overlap  (no turns)
simulations(2).run(9).targets = [1 8];

simulations(3).environment = OPT_TESSELATION_SQUARE_LINE;
simulations(3).maxRounds = 2500;
start = 0;
for tmp = start + 1 : 16
    simulations(3).run(tmp - start).sources = [1];
    simulations(3).run(tmp - start).targets = [tmp];
end

simulations(4).environment = OPT_TESSELATION_EQTRI_LINE;
simulations(4).maxRounds = 2500;
start = 0;
for tmp = start + 1 : 16
    simulations(4).run(tmp - start).sources = [1];
    simulations(4).run(tmp - start).targets = [tmp];
end

simulations(5).environment = OPT_TESSELATION_RTRI_LINE;
simulations(5).maxRounds = 2500;
start = 0;
for tmp = start + 1 : 16
    simulations(5).run(tmp - start).sources = [1];
    simulations(5).run(tmp - start).targets = [tmp];
end



simulations(7).environment = OPT_TESSELATION_SQUARE;
simulations(7).maxRounds = 2500;
simulations(7).run(1).sources = [3 18]; % square: overlapping
simulations(7).run(1).targets = [23 8];

simulations(8).environment = OPT_TESSELATION_SQUARE;
simulations(8).maxRounds = 2500;
simulations(8).run(1).sources = [3 11]; % square left-right throughput, 2 overlap (has turns)
simulations(8).run(1).targets = [23 19]; 
simulations(8).run(1).failed = [14 9 7];

simulations(9).environment = OPT_TESSELATION_SQUARE;
simulations(9).maxRounds = 2500;
simulations(9).run(1).sources = [3 7]; % square left-right throughput, 3 overlap  (has turns)
simulations(9).run(1).targets = [23 19]; 
simulations(9).run(1).failed = [12 17 9 14];


simulations(10).environment = OPT_TESSELATION_SQUARE;
simulations(10).maxRounds = 2500;
simulations(10).run(1).sources = [3 18 25 21 6]; % square: overlapping, multi-lock disjoint case
simulations(10).run(1).targets = [23 8 5 1 16];

simulations(11).environment = OPT_TESSELATION_SQUARE;
simulations(11).maxRounds = 2500;
simulations(11).run(1).sources = [4 20 22 6]; % square corners
simulations(11).run(1).targets = [25 21 1 5];
simulations(11).run(1).failed = [7 8 9 12 13 14 17 18 19];

simulations(11).run(2).sources = [4 20 22 6 25 21 1 5]; % square corners, 2 copies
simulations(11).run(2).targets = [25 21 1 5 4 20 22 6];
simulations(11).run(3).failed = [7 8 9 12 13 14 17 18 19];

simulations(11).run(3).sources = [4 20 22 6 4 20 22 6 4 20 22 6]; % square corners, 3 copies
simulations(11).run(3).targets = [25 21 1 5 25 21 1 5 25 21 1 5];
simulations(11).run(3).failed = [7 8 9 12 13 14 17 18 19];


simulations(12).environment = OPT_TESSELATION_SQUARE;
simulations(12).maxRounds = 2500;
simulations(12).run(1).sources = [10 16 24 2]; % square: normal intersection, all straight through
simulations(12).run(1).targets = [6 20 4 22];
simulations(12).run(1).failed = [3 13 23];

simulations(12).run(2).sources = [16 24 10 2]; % square: normal intersection, all right turns
simulations(12).run(2).targets = [22 20 4  6];
simulations(12).run(2).failed = [3 13 23 21 25 5 1];

simulations(12).run(3).sources = [16 24 10 2]; % square: normal intersection, all left turns
simulations(12).run(3).targets = [4  6  22 20];
simulations(12).run(3).failed = [3 13 23 21 25 5 1];


simulations(13).environment = OPT_TESSELATION_SNUB_SQUARE_TILING;
simulations(13).maxRounds = 250;
simulations(13).run(1).sources = [3 4 6 5]; % snub square, 4 copies
simulations(13).run(1).targets = [10 15 24 23];% snub square, 4 copies
%targets = [10 15 11 18];% snub square, 3 copies
%sources = [3 4 6 5]; % snub square, 3 copies


simulations(14).environment = OPT_TESSELATION_PARALLELOGRAM;
simulations(14).maxRounds = 2500;
simulations(14).run(1).sources = 4; % parallelogram
simulations(14).run(1).targets = 1; % parallelogram with failures
simulations(14).run(1).failed = [3 7 11 13 9 5]; % parallelogram

simulations(15).environment = OPT_TESSELATION_PARALLELOGRAM;
simulations(15).maxRounds = 2500;
simulations(15).run(1).sources = [3 11 7 19]; % parallelogram
simulations(15).run(1).targets = [23 15 1 25]; % parallelogram with failures


simulations(16).environment = OPT_TESSELATION_SQUARE;
simulations(16).maxRounds = 2500;
simulations(16).run(1).sources = [2 24]; % one-lane bridge
simulations(16).run(1).targets = [22 4]; % one-lane bridge, squares
simulations(16).run(1).failed = [9 14 19 7 12 17   1  6 11 16 21 5 10 15 20 25]; % one-lane bridge

simulations(16).run(2).sources = [2 24]; % one-lane bridge
simulations(16).run(2).targets = [22 4]; % one-lane bridge, squares
simulations(16).run(2).failed = [9 14 19]; % one-lane bridge, dynamic failure
simulations(16).run(2).faildynamic = [7 12 17   1  6 11 16 21 5 10 15 20 25];


simulations(17).environment = OPT_TESSELATION_SQUARE;
simulations(17).maxRounds = 200;
simulations(17).run(1).sources = [5 11 21]; % source fairness assumptions
simulations(17).run(1).targets = [25 12 1];
simulations(17).run(1).failed = [];
simulations(17).run(1).faildynamic = [20 19 14 9];


simulations(18).environment = OPT_TESSELATION_SQUARE;
simulations(18).maxRounds = 500;
simulations(18).run(1).sources = [10 16]; % one-lane bridge (vertical)
simulations(18).run(1).targets = [6 20]; % one-lane bridge, (vertical)
simulations(18).run(1).failed = [1 2 3 4 5 7 8 9 17 18 19 21 22 23 24 25]; % one-lane bridge (vertical)


simulations(19).environment = OPT_TESSELATION_SQUARE;
simulations(19).maxRounds = 200;
simulations(19).run(1).sources = [2  8 11 20 24 5]; % lockcolors example
simulations(19).run(1).targets = [17 6 13 18 14 15]; % lockcolors example
simulations(19).run(1).failed = [];
simulations(19).run(1).faildynamic = [];


%targets = 10;  % hexagon
%sources = 16; % hexagon

%sources = [3 11]; % square left-right
%targets = [23 15]; 






%targets = [25 21 1 5 23 15]; % square corners with more
%sources = [4 20 22 6 2 16];
%failed = [7 8 9 12 13 14 17 18 19];







if ~opt_single
    min_sim = 1;
    max_sim = length(simulations);
else
    min_sim = opt_single;
    max_sim = opt_single;
    min_run = opt_run;
    max_run = opt_run;
end

for iSim = min_sim : max_sim
    
    if ~opt_run
        min_run = 1;
        max_run = length(simulations(iSim).run);
    end
    
    for iRun = min_run : max_run
        
        sources = simulations(iSim).run(iRun).sources;
        targets = simulations(iSim).run(iRun).targets;
        opt_tesselation_type = simulations(iSim).environment;
        try
            failed = simulations(iSim).run(iRun).failed;
        catch
            failed = [];
        end
        
        try % need to be separate blocks, each can fail independently
            faildynamic = simulations(iSim).run(iRun).faildynamic;
        catch
            faildynamic = [];
        end
        
        opt_triangulation = 0; % 0 = don't use delauny triangulation (use a uniformly generated set of vertices instead, for throughput simulations)
        
        placedColors = [];
        lockedColors = [];

epsPlot = 0.25;
epsIllegal = 0.0001;

L = 0.1;
%L = 0.05;       % entity radius; for throughput
rs = L / 10;    % safety radius
d = L + rs; 
m = 2;          % number of dimension (use this variable when necessary to refer to the dimensions)
F = 0;          % number of cells to fail
safeDistance = 3 * L + rs;

NT = 1; % number of targets
NS = NT; % number of sources

MAX_ROUNDS = simulations(iSim).maxRounds;

% options
mode_debug = 0;
opt_randic = 0; % random initial condition?
opt_fixedic = 1; % case number for fixed initial conditions (see switch/case statement below); only applies if opt_randic = 0
opt_center_entity = 0; % 1 = place one entity at center of cell at a time, 0 = place (about) numToTryPlacing
%numToTryPlacing = ceil(1/L)/2; % todo: generalize and find side length
numToTryPlacing = 20;
opt_roundsToPlace = 1; % try to place entities every opt_roundsToPlace rounds
opt_display = 1;        % don't draw system evolution if 0, draw if 1
%opt_display_update = 10;
opt_display_update = 1;
opt_display_progress = 100;  % number of rounds to display a text update on progress if we aren't drawing the system evolution
opt_decimal_tolerance = -4; % decimal tolerance (10^opt_decimal_tolerance)
opt_badCell = 1;
opt_singleEntity = 0;
opt_video = 0; % record video

if opt_video
    writerObj = VideoWriter('factory.avi');
    writerObj.FrameRate = 5;
    open(writerObj);
end

opt_linprog_options=optimset('Display', 'off'); % used to supress messages from linear programming solver

opt_routing_disjoint = 0; % choose routing mode: we may in the future want to implement the disjoint paths algorithm, but for now will not

gcf;
%set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
axis equal;
set(gcf, 'Position', [1 1 800 800]); % Maximize figure.

set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer');

% set failed cells manually
failround = 75;
%failround = -1;
opt_faildynamic = 1; % option to allow failures occurring while running


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
        case OPT_TESSELATION_RECTANGULAR
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

% must be defined after diam
opt_firstFullRound = 2*diam; % first full round to place entities, etc.; in general needs to be 2*N

if opt_triangulation
    T = (1:N)';
    neigh = tr.neighbors(); % set of all neighbors
    cc = tr.circumcenters(); % set of all circumcenters
    [ic,icr] = tr.incenters(); % set of all incenters with corresponding radii (inscribed circle center and radius)
    idx1 = T < neigh(:,1); % modify indexing to make appropriate neighbor graph
    idx2 = T < neigh(:,2);
    idx3 = T < neigh(:,3);
    neigh = [T(idx1) neigh(idx1,1); T(idx2) neigh(idx2,2); T(idx3) neigh(idx3,3)]';
end

% generate source and target pairs
% identifiers of the source cells: note, must be same length as targets, as each source will produce entities of the same type of the corresponding entry of the targets vector
if opt_randic
    failed = ceil(N*rand(F,1));
    targets = ceil(N*rand(NT,1));
    sources = ceil(N*rand(NS,1));
    
    % ensure targets of colors are unique (not required for analysis, just more interesting simulations: if the targets are the same cell, then just about everything follows from the original ICDCS paper)
    while length(unique(targets)) ~= NT || ~isempty(intersect(targets,failed))
       targets = ceil(N*rand(NT,1));
    end
    
    % ensure targets, failed, and sources are disjoint and sources are also unique (again, not required for analysis, just for more interesting simulations)
    while ~isempty(intersect(targets,sources)) || ~isempty(intersect(sources,failed)) || length(unique(sources)) ~= NS
       sources = ceil(N*rand(NS,1));
    end
else
    %failed = [];
    % identifiers of the sinks (target cells)
    switch opt_fixedic
        case 1
            targets = targets;
            sources = sources;
        case 2
            targets = targets;
            sources = sources;
            
        % test cases on uniform grid
        case 3
            targets = [32; 3]; % nice test
            sources = [23; 6]; % nice test
        case 4
            targets = [32; 23];
            sources = [3; 6];
        case 5
            targets = [32; 3; 26; 17]; % nice test, 4 colors: some deadlocks
            sources = [23; 6; 2; 10]; % nice test, 4 colors: some deadlocks
        case 6
            targets = [32; 3; 9; 20]; % 2 disjoint locksets
            sources = [23; 6; 12; 31]; % 2 disjoint locksets
        case 7
            targets = [22; 16; 4; 14; 3; 6]; % single large lockset
            sources = [19; 5; 13; 10; 23; 32]; % one large lockset
        case 8
            targets = [32; 3; 9; 20; 14; 16]; % 3 disjoint locksets
            sources = [23; 6; 12; 31; 13; 19]; % 3 disjoint locksets
        case 9
            targets = [18 5];
            sources = [1 11];
    end
end

% fix NT if we manually set fixed targets
if NT ~= length(targets)
    NT = length(targets);
end

% fix NS if we manually set fixed sources
if NS ~= length(sources)
    NS = length(sources);
end

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

for ti = 1 : NT
    path(ti).lockcolors = [];
    path(ti).int = [];
    path(ti).intEmpty = [];
    path(ti).intColor = [];
    path(ti).donecolors = [];
    path(ti).nlock = 0;
end

% fail F cells randomly
% for i = 1 : F
%     f = ceil(rand(1,1)*N);
% 
%     % don't fail targets or sources
%     if intersect([targets; sources], f)
%         continue;
%     end
%     failed = [failed; f];
% end

% set distance and next pointers for failures
if ~isempty(failed)
    for i = 1 : length(failed)
        f = failed(i);
        Cell(f).failed = 1;
        for tt = 1 : NT
            Cell(f).path(tt).prev = [];
            Cell(f).next(tt) = f;
            Cell(f).dist(tt) = inf;
        end
    end
end

% error checking
if sum(targets > N) ~= 0 || sum(sources > N) ~= 0 || sum(failed > N) ~= 0 || sum(faildynamic > N) ~= 0
    'ERROR: bad source, target, or failed IDs entered (greater than N)'
    return
end

% find neighbors (inefficient all to all compare, could be done smarter)
for i = 1 : N
    for j = 1 : N
        % cells are neighbors iff they share an edge
        if i ~= j && size(intersect(roundn(Cell(i).vertices,opt_decimal_tolerance), roundn(Cell(j).vertices,opt_decimal_tolerance),'rows'),1) >= 2
            Cell(i).nbrs = unique([Cell(i).nbrs, j]);
            %Cell(j).nbrs = [Cell(j).nbrs, i];
        end
    end
    Cell(i).nfnbrs = Cell(i).nbrs;
end

% associate neighbors with sides
% needs to be a separate loop: relies on Cell(j).vertices being constructed for j \in nbrs(i)
for i = 1 : N
    for s = 1 : Cell(i).num_sides
        Cell(i).side(s).safe = 1;   % side s starts safe (we run safety check before moving)
        
        for j = Cell(i).nbrs
            % delete bad neighbor indices
            if j > N || j <= 0
                Cell(i).nbrs = Cell(i).nbrs(Cell(i).nbrs ~= j); % remove neighbors with too large or too small indices
                continue;
            elseif isnan(j)
                Cell(i).nbrs = Cell(i).nbrs(~isnan(Cell(i).nbrs)); % remove NaN neighbors
                continue;
            end
            
            if size(intersect(roundn(Cell(i).side(s).vertices, opt_decimal_tolerance), roundn(Cell(j).vertices, opt_decimal_tolerance), 'rows')) == [2 2] % share two vertices => share an edge => adjacent triangles and this side is not on the boundary of the triangulation
                Cell(i).side(s).nbr = j;    % index of neighbor which touches cell i on side s
                Cell(i).side(s).boundary = 0;
            end
        end
    end
end

% Initialize all source cells
for i = 1 : length(sources)
    Cell(sources(i)).color = i;                             % set the "color" / type to be equal to the entry in the targets vector, i.e., 1, 2, 3, etc.
end

% Initialize all target cells to be 0 distance away from the target
for i = 1 : length(targets)
    Cell(targets(i)).dist(i) = 0;   % set the distance to 0 away from this color
    Cell(targets(i)).color = i;     % set the "color" / type to be equal to the entry in the targets vector, i.e., 1, 2, 3, etc.
end

opt_waitplot = 1; % on = 1, off = 0
waitplot = 0.000001; % wait (in seconds) between plots

indexEntity = 1; % global index for entities (to see if we duplicate, compare them between rounds, etc.)

CellOld = Cell; % copy the current cell information (models communication being delayed a round)
                % To model communications, from cell i, any access to
                % information of neighbors must be from the CellOld
                % variable, which represents what was communicated to i
                % from a neighbor at the last round.  For the current
                % round, update the Cell variable.
locks = [];
overlap = 0;
startOver = ones(NT,1);

% main simulation loop, k is the round index
for k = 2 : MAX_ROUNDS
    NF = [1:N];
    for i = 1 : N
        if Cell(i).failed
            NF = setdiff(NF, i);
        end
        
        Cell(i).nfnbrs = Cell(i).nbrs;
        % could use NF, but then would need another loop
        for j = Cell(i).nbrs
            if Cell(j).failed
                Cell(i).nfnbrs = setdiff(Cell(i).nfnbrs, j);
            end
        end
    end
    
    for ti = 1 : NT
        if isempty(path(ti).int)
            placedColors = setdiff(placedColors, ti);
        end
    end
    
    for i = NF
        % compute prev
        for ti = 1 : NT
            for j = Cell(i).nfnbrs
                Cell(i).path(ti).prev = [];
                if Cell(j).next(ti) == i
                    Cell(i).path(ti).prev = unique([Cell(i).path(ti).prev; j]);
                end
            end
        end
    end
    
    % fail a cell after we've started running (for plots)
    if k == failround && opt_faildynamic
        failed = [failed, faildynamic];
    
        % set distance and next pointers for failures
        if ~isempty(failed)
            for i = 1 : length(failed)
                f = failed(i);
                Cell(f).failed = 1;
                for tt = 1 : NT
                    Cell(f).prev(tt) = f;
                    Cell(f).next(tt) = f;
                    Cell(f).dist(tt) = inf;
                end
            end
        end
    end

%     for ti = 1 : NT
%         for i = NF
%             if sources(ti) == i && ~isempty(intersect(placedColors, ti)) % this color already put stuff on the source, block it
%                 %Cell(i).path(ti).rlock = 1;
%             else
%                 %Cell(i).path(ti).rlock = isempty(Cell(i).Entities);
%             end
%         end
%     end
%     
%     for ti = 1 : NT
%         for i = NF
%             if ~(sources(ti) == i && ~isempty(intersect(placedColors, ti))) % this color already put stuff on the source, block it
%                 for j = Cell(i).nfnbrs
%                     %Cell(i).path(ti).rlock = Cell(i).path(ti).rlock & Cell(j).path(ti).rlock;
%                 end
%             end
%         end
%     end


    for ti = 1 : NT
        path(ti).intEmpty = 1; % if all paths empty, then pick a color
        path(ti).nonempty = 0;
        path(ti).done = 0;
    end
    
    % need separate loop to get lockcolors correct (uses path.int)
    for i = NF
        for ti = 1 : NT
            if ~isempty(intersect(Cell(i).path(ti).int, i))
                % determine if lockset (intersection) empty or not
                path(ti).intEmpty = path(ti).intEmpty & isempty(Cell(i).Entities); % note: not indexed by cell, global
                if ~isempty(Cell(i).Entities)
                    path(ti).intColor = ti;
                end
                path(ti).nlock = 1;
            end
            
            if ~isempty(intersect(Cell(i).path(ti).ids, i))
                path(ti).nonempty = path(ti).nonempty | ~isempty(Cell(i).Entities); % note: not indexed by cell, global
            end
        end
    end
    
    reset = ones(NT,1);
    for ti = 1 : NT
        for i = NF
            % only do this for cells on the path intersection
            if ~isempty(intersect(Cell(i).path(ti).ids, i))                
                reset(ti) = reset(ti) & path(ti).intEmpty  & ~isempty(intersect(placedColors, ti)) & ~isempty(intersect(lockedColors, ti)) & ~startOver(ti); % compute reset as the conjunction of all rlock variables for each color
            end
            %& ~(Cell(i).detectNextFailed > 0)
            
            % failure detection (to avoid deadlocking after failures)
            if (Cell(i).detectNextFailed > 0 && k < Cell(i).detectNextFailed + opt_firstFullRound)
                %locks = [];
                
                if ~path(ti).intEmpty
                    locks = setdiff(locks, path(ti).lockcolors);
                    locks = unique([locks; path(ti).intColor]);
                end
                %break;
            end
        end
    end
    
    % mutual exclusion between overlapping colors
    if k > opt_firstFullRound || length(targets) == 1 % must wait at least 2*N rounds for information to definitely have propogated between every cell (most cases would allow about 2*sqrt(N), but 2*N is the case of a line graph)
        for ti = 1 : NT
            if ~path(ti).done && path(ti).nlock
            
                if isempty(intersect(path(ti).lockcolors, ti))
                    placedColors = setdiff(placedColors, ti);
                end

                for i = NF
                    % change the lock
                    if (path(ti).intEmpty) || ( ~isempty(Cell(i).Entities) && startOver(ti) ) || (~isempty(Cell(i).path(ti).nlock) && (startOver(ti) || ~isempty(intersect(sources, i)) || ~isempty(Cell(i).Entities) )) % cell i has entities or is a source
                        if ( allLockIntEmpty(locks, path) || reset(ti) || isempty(intersect(locks, Cell(i).path(ti).lockcolors))) && path(ti).intEmpty % no one has lock, or reset this lock
                            % this is the color that has the lock and is losing it
                            rc = Cell(i).path(ti).lockcolors( find( reset(Cell(i).path(ti).lockcolors) == 1) );
                            
                            if ~isempty(intersect(rc, ti))                  
                                locks = setdiff(locks, rc); % ensure take back lock

                                % all locks are reset, start over
                                if ~isempty(intersect(placedColors, rc)) && length( rc ) == length( Cell(i).path(ti).lockcolors) && isempty( find( rc == Cell(i).path(ti).lockcolors == 0) )
                                    startOver(ti) = 1;
                                    lockedColors= setdiff(lockedColors, rc); % todo: rc?
                                    placedColors = setdiff(placedColors, rc);
                                end
                            else
                                if mode_debug
                                    locks
                                end
                                locks = setdiff(locks, path(ti).lockcolors); % ensure no lock colors in set of locks (otherwise violate mutex)
                                locks = unique([locks; ti]);

                                if startOver(ti)
                                    startOver(ti) = 0;
                                end
                                
                                lockedColors = unique([lockedColors; ti]);
                                
                                path(ti).done = 1;
                                for tj = path(ti).lockcolors'
                                    path(tj).done = 1;                                    
                                end

                                ['Lock changed: added ', num2str(ti), ' at round ', num2str(k)]
                                if mode_debug
                                    locks
                                    reset
                                    rc
                                    placedColors
                                    startOver
                                    path.done
                                end

                                break; % break out of i loop
                                %pause
                            end
                        end
                    end
                end
            end
        end
    end

    % every X rounds, try to place lots of entities on sources safely
    % don't place in first round because that could cause deadlocks
    if mod(k,opt_roundsToPlace) == 0 && (k >= opt_firstFullRound || length(targets) == 1)
        for tt = 1 : NS
            s = sources(tt); % source cell id
            j = 0; % index of number of entities to try placing
            safe = 1;

            % place a single entity at the centroid (for throughput measurements)
            if opt_center_entity || Cell(s).safeRegion.R < 0
                x = Cell(s).centroid;
                for p = 1 : length(Cell(s).Entities)
                    % unsafe if either too close, or trying to put this color entities on a source cell which already has another color on it (can occur, e.g., red on green source, then green source cannot place new entities)
                    if norm( Cell(s).Entities(p).x - x, 2) < 2*d || Cell(s).Entities(p).color ~= tt
                        safe = 0;
                        break; % short-circuit to exit loop since safety is for all entities
                    end
                end
                % only add entities to sources if it's safe to do so
                % * Don't need lock, or need the lock and have it
                if safe && isempty(intersect(placedColors, tt)) && (isempty(Cell(s).path(tt).nlock) || (~isempty(Cell(s).path(tt).nlock) && Cell(s).path(tt).lock && ~isempty(intersect(locks, tt)) && path(tt).intEmpty))
                    % todo: check, removed
                    en = length(Cell(s).Entities);
                    Cell(s).Entities(en + 1).x = x;
                    Cell(s).Entities(en + 1).id = indexEntity;
                    Cell(s).Entities(en + 1).color = tt;             % set the color of the entity to be the color of the source
                    Cell(s).Entities(en + 1).moved = 1;             % set a bit that this has already moved this round
                    indexEntity = indexEntity + 1;                  % increment global index since we've created an entity
                end
                
                if (~isempty(Cell(s).Entities) || length(Cell(s).Entities) > 0) && Cell(s).Entities(1).color == tt
                    placedColors = unique([placedColors; tt]);
                end
            else % try to place a bunch of entities safely each round
                while j <= numToTryPlacing
                    safe = 1; % todo: why was this commented?
                    % generate a random point in the source (convex combination of the vertices => always inside the convex set)
                    % note: we do this inside the smaller homothetic triangle to ensure the entities don't start too close to a face (i.e., where they would need to transfer already)
                    ccSides = rand(Cell(s).num_sides,1); % need a point from each side
                    while sum(ccSides(1: Cell(s).num_sides - 1)) >= 1
                        ccSides = rand(Cell(s).num_sides,1); % need a point from each side
                    end
                    ccSides( Cell(s).num_sides ) = 1 - sum(ccSides(1: Cell(s).num_sides - 1)); % ensures cc1 + cc2 + cc3 + ... = 1 and cc1,cc2,cc3,... > 0 (convex combination)
                    x = [0, 0];
                    for side = 1 : Cell(s).num_sides
                        x = x + ccSides(side,:) * Cell(s).illegalRegion.vertices(side,:);
                    end

                    for p = 1 : length(Cell(s).Entities)
                        % unsafe if either too close, or trying to put this color entities on a source cell which already has another color on it (can occur, e.g., red on green source, then green source cannot place new entities)
                        if norm( Cell(s).Entities(p).x - x, 2) <= 2.5*d || Cell(s).Entities(p).color ~= tt
                            safe = 0;
                            break; % short-circuit to exit loop since safety is for all entities
                        end
                    end
                    % only add entities to sources if it's safe to do so
                    % * Don't need lock, or need the lock and have it
                    if safe && isempty(intersect(placedColors, tt)) && (isempty(Cell(s).path(tt).nlock) || (~isempty(Cell(s).path(tt).nlock) && Cell(s).path(tt).lock && ~isempty(intersect(locks, tt)) && path(tt).intEmpty))
                        %~Cell(s).path(tt).rlock)
                        %&& ~isempty(intersect(locks, tt)))) 
                        % todo: check, removed && Cell(s).path(i).lock)
                        en = length(Cell(s).Entities);
                        Cell(s).Entities(en + 1).x = x;
                        Cell(s).Entities(en + 1).id = indexEntity;
                        Cell(s).Entities(en + 1).color = tt;             % set the color of the entity to be the color of the source
                        Cell(s).Entities(en + 1).moved = 1;             % set a bit that this has already moved this round
                        indexEntity = indexEntity + 1;                  % increment global index since we've created an entity
                        
                        %Cell(s).path(tt).rlock = 1; % fairness: reset lock for source after placing
                    end
                    j = j + 1;
                end

                if (~isempty(Cell(s).Entities) || length(Cell(s).Entities) > 0) && Cell(s).Entities(1).color == tt
                    placedColors = unique([placedColors; tt]);
                end
            end
        end
    end

    % color cells based on entity types
    for i = NF
        Cell(i).etype = [];
        for p = 1 : length(Cell(i).Entities)
            if isempty(intersect(Cell(i).etype, Cell(i).Entities(p).color))
                Cell(i).etype = [Cell(i).etype Cell(i).Entities(p).color];
            end
        end
        Cell(i).etype = unique(Cell(i).etype); % ensure unique
        
        % if a cell doesn't have a previous entity type yet, go ahead and update
        % otherwise, it is only ever updated by the move function when the last entity leaves the cell
        if isempty(Cell(i).lastEtype)
            Cell(i).lastEtype = Cell(i).etype;
        end
        
        for tt = 1 : NT
            % see if next pointer is failed, detect failure
            np = Cell(i).next(tt);
            if isinf(Cell( np ).dist(tt)) && k >= opt_firstFullRound
                Cell(i).detectNextFailed = k;
            end
        end
    end
    
    % iterate over all the targets (create a different route for different color types; note that we are assuming there is a unique "color" or type corresponding to each target)
    for tt = 1 : NT
        % iterate over all the non-faulty cells
        for i = NF
            Cell(i).path(tt).int = []; % reset to start empty
            
            % don't change sinks' next or dist for their own color (but do for other colors)
            if targets(tt) == i
                continue;
            end

            % ROUTING
            % iterate over Cell i's neighbors
            % reset these sets
            shortest = inf;
            shortestid = inf;

            for j = CellOld(i).nfnbrs
                % disjoint path version: this won't stabilize: adds an extra cost for going over other edges, but this graph is a spanning tree, thus all edges will overlap
%                 distjPath = 0;
%                 for tj = NT : NT
%                     if tj ~= tt && Cell(i).next(tj) ~= i && Cell(i).next(tt) ~= i && ...
%                             (~isempty(intersect(Cell(i).path(tj).prev, Cell(i).next(tt)))) && ...
%                             (~isempty(intersect(Cell(i).path(tt).prev, Cell(i).next(tj))))
%                         %(isempty(Cell(i).path(tt).prev) || ~isempty(intersect(Cell(i).path(tt).prev, Cell(i).next(tj))))
%                         %(isempty(Cell(i).path(tj).prev) || ~isempty(intersect(Cell(i).path(tj).prev, Cell(i).next(tt)))) && ...
%                         distjPath = 1
%                     end
%                 end
                %  + distjPath
                if CellOld(j).dist(tt) < min(shortest, CellOld(i).dist(tt)) && j < shortestid % basic search
                    % todo: use argmin over ids? should prevent this non-contiguous lockset we saw
                    shortest = CellOld(j).dist(tt); % todo / note: gets different result using CellOld, needs concurrent computation; will be same length, but different ids
                    shortestid = j;
                end
            end

            Cell(i).dist(tt) = shortest + 1;

            if isinf(shortestid)
                Cell(i).next(tt) = i; % point at self if didn't find a shorter one
            else % avoid already tried ones to get disjoint paths
                Cell(i).next(tt) = shortestid;
            end
        end
    end
    
    % compute the path
    for i = NF
        for tt = 1 : NT
            if opt_routing_disjoint == 0
                % inductively defined as:
                % path_0: if i has entities or is a source
                % path_{i}: if i is in the ids of this color, add its' next

                % ensure routes have stabilized (sort of...)
                %if CellOld(i).next ~= Cell(i).next
                %    break; % don't do anyone if any have changed
                %end

                if ((~isempty(Cell(i).Entities) && Cell(i).etype == tt) || sources(tt) == i || ~isempty(intersect(Cell(i).path(tt).ids, i))) && Cell(i).next(tt) ~= i
                    Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; i; Cell(i).next(tt)]);
                end

                % todo: after a long time from the last failure, we may be able to remove cells from
                % being on the path, e.g., after there are only entities on
                % cells between the source and path, we can remove the old ones

                % aggregate paths of all neighbors (gossip)
                for j = Cell(i).nbrs
                    Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; CellOld(j).path(tt).ids]); % note use of old
                end
            else
                restart = 1;
                % reversal version
                while (restart)
                    restart = 0;
                    if ((~isempty(Cell(i).Entities) && Cell(i).etype == tt) || sources(tt) == i || ~isempty(intersect(Cell(i).path(tt).ids, i))) && Cell(i).next(tt) ~= i
                        Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; i; Cell(i).next(tt)]);
                    end

                    for tj = 1 : NT
                        if tt ~= tj
                            % path for tj subseteq tt
                            common = intersect(Cell(i).path(tt).ids,Cell(i).path(tj).ids);
                            if ~isempty(common) && length(common) == length(Cell(i).path(tj).ids) && isempty(find((common == Cell(i).path(tj).ids) == 0))
                                %'subset'
                                %for rev = 1 : length(Cell(i).path(tt).ids)
                                %    ri = Cell(i).path(tt).ids(rev);
                                %    tmp = ceil(length(Cell(ri).nfnbrs)*rand(1,1));
                                %    Cell(ri).next(tt) = Cell(ri).nfnbrs(tmp);
                                %end
                            end
                        end
                    end

                    % todo: after a long time from the last failure, we may be able to remove cells from
                    % being on the path, e.g., after there are only entities on
                    % cells between the source and path, we can remove the old ones

                end
                % aggregate paths of all neighbors (gossip)
                for j = Cell(i).nbrs
                    Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; CellOld(j).path(tt).ids]); % note use of old
                end
            end

%             % a-star version
%             restart = 1;
%             retry = 0;
%             nextToTry = Cell(i).next(tt);
%             while restart
%                 %Cell(i).path(tt).ids = CellOld(i).path(tt).ids;
%                 restart = 0;
%                 if ((~isempty(Cell(i).Entities) && Cell(i).etype == tt) || sources(tt) == i || ~isempty(intersect(Cell(i).path(tt).ids, i))) && nextToTry ~= i && isempty(intersect(Cell(i).path(tt).nextTried, nextToTry))
%                     if retry
%                         Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; nextToTry]);
%                     else
%                         Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; i; nextToTry]);
%                     end
%                 end
% 
%                 % todo: after a long time from the last failure, we may be able to remove cells from
%                 % being on the path, e.g., after there are only entities on
%                 % cells between the source and path, we can remove the old ones
%                  
%                 Cell(i).path(tt).int = []; % reset to start empty
%                  
%                 for tj = 1 : NT
%                     if tt == tj
%                        continue
%                     end
%                     pint = unique(intersect( Cell(i).path(tt).ids, Cell(i).path(tj).ids ));
%                     Cell(i).path(tt).int = unique([Cell(i).path(tt).int; pint]);
%                 end
%  
%                 if ~isempty(Cell(i).path(tt).int) && ~isempty(intersect(Cell(i).path(tt).int, i))
%                     'retry'
%                     i
%                     tt
%                     
%                     Cell(i).path(tt).nextTried = unique([Cell(i).path(tt).nextTried; nextToTry]);
%                     
%                     tmpnbrs = Cell(i).nfnbrs;
% %                     tmpnbrs = [];
% %                     for j = Cell(i).nfnbrs
% %                         %if Cell(j).next(tt) ~= i
% %                         if ~isempty(intersect(Cell(i).path(tt).ids, j))
% %                             tmpnbrs = [tmpnbrs; j];
% %                         end
% %                     end
% 
%                     tmp = setdiff(tmpnbrs, Cell(i).path(tt).nextTried);
%                     if isempty(tmp) % have to quit
%                         restart = 0;
%                         if i == 17 || i == 12 || i == 22 || i == 7 || i == 2
%                             Cell(i).path(tt).nextTried
%                             %pause
%                         end
%                     else
%                         len = length(tmp);
%                         Cell(i).path(tt).ids = setdiff(Cell(i).path(tt).ids, [nextToTry]);
%                         %Cell(i).path(tt).ids = [i];
%                         %nextToTry = tmp(1);
%                         nextToTry = tmp( ceil(len*rand(1,1)) );
%                         
%                         if i == 17 || i == 12 || i == 22 || i == 7 || i == 2
%                             'trying next'
%                             nextToTry
%                             %pause
%                         end
% 
%                         restart = 1;
%                         retry = 1;
%                     end
%                 end
%             end
%             
%             % aggregate paths of all neighbors (gossip)
%             for j = Cell(i).nfnbrs
%                 
%                 % if i isn't in it's own path, we need to not gossip it around
%                 %if isempty(intersect(Cell(i).path(tt).ids, i)) && ~isempty(intersect(CellOld(j).path(tt).ids, i))
%                 %    Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; CellOld(j).path(tt).ids]); % note use of old
%                 %    Cell(i).path(tt).ids = setdiff(Cell(i).path(tt).ids, i);
%                 %else
%                 %    Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; CellOld(j).path(tt).ids]); % note use of old
%                 %end
%                 Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; CellOld(j).path(tt).ids]); % note use of old
%             end
        end
    end
    
    % compute the lockset
    for i = NF
        for ti = 1 : NT
            Cell(i).path(ti).int = []; % reset to start empty
            for tj = 1 : NT
                if ti == tj
                    continue
                end
                
                % all pairs intersection
                pint = unique(intersect( Cell(i).path(ti).ids, Cell(i).path(tj).ids ));
                
                % if there exists a single color in pint with a different next pointer in pint, then all colors on that cell must be locked
                % * This allows us to lock "directions" as opposed to colors (or to lock colors, but allow any color with that direction to move)
                %if ~isempty(pint) % && Cell(i).next(ti) ~= Cell(i).next(tj) % can't do the same nexts always, would improve throughput, but have counterexample
                    Cell(i).path(ti).int = unique([Cell(i).path(ti).int; pint]);
                    % if i is in the intersection, and the two colors do not point in the same direction, then it will need to be locked
                    if ~isempty(intersect(Cell(i).path(ti).int, i))  %todo: was pint
                        % Cell i needs the locks for colors ti and tj
                        Cell(i).path(ti).nlock = unique([Cell(i).path(ti).nlock; ti; tj]);
                    end
                %end
            end
        end
    end
    
    % need separate loop to get lockcolors correct (uses path.int)
    for i = NF
        for ti = 1 : NT
            for tj = 1 : NT
                %if ti ~= tj %&& ~isempty(intersect(Cell(i).path(ti).int, Cell(i).path(tj).int))
                %if ti ~= tj
                if ~isempty(intersect(Cell(i).path(ti).int, Cell(i).path(tj).int))                    
                    % union of all overlapping colors define the "lockset"
                    Cell(i).path(ti).lockset = unique([Cell(i).path(ti).lockset; Cell(i).path(ti).int; Cell(i).path(tj).int]);
                    Cell(i).path(ti).lockcolors = unique([Cell(i).path(ti).lockcolors; Cell(i).path(tj).lockcolors; ti; tj]); % note, if we don't add also tj.lockcolors, we can actually allow higher throughput, as this would be the case where say 4 paths intersect at some points, but only 3 at all the same points simultaneuosly (but this may lead to deadlocks depending on how we cycle through the colors)
                    path(ti).int = Cell(i).path(ti).int;
                    path(ti).lockset = Cell(i).path(ti).lockset;
                    path(ti).lockcolors = Cell(i).path(ti).lockcolors;
                end
            end
        end
    end

    % LOCKING
    for i = NF
        for ti = 1 : NT
            % cell i needs a lock for this path
            if ~isempty(intersect(Cell(i).path(ti).ids, i)) && ~isempty(intersect(Cell(i).path(ti).nlock, ti)) && ~isempty(intersect(locks, ti))
                if ~Cell(i).path(ti).rlock || ~isempty(Cell(i).Entities) % don't reset lock
                    Cell(i).path(ti).lock = 1;
                else
                    %['lock really reset for cell ', num2str(i)]
                    %Cell(i).path(ti).lock = 0;
                end
            else
                Cell(i).path(ti).lock = 0;
            end
        end
    end
    
    % reset safety and signals
    for i = NF
        % start all signals off
        Cell(i).signal = i;
        for sd = 1 : Cell(i).num_sides
            Cell(i).side(sd).safe = 1; % start all sides safe
        end
    end
    
    for i = NF
        %checking if the current cell has any safe sides to transfer
        if ~isempty(Cell(i).Entities)
            for j = Cell(i).nbrs
                % find the side corresponding to neighbor j (otherwise we don't have to do anything, since a side without a neighbor never has entities coming from it, nor going to it)
                sd = [];
                for stmp = 1 : Cell(i).num_sides
                    if j == Cell(i).side(stmp).nbr
                        sd = stmp;
                        break;
                    end
                end
                
                if isempty(sd)
                   continue;
                end
                
                % see if there are any entities in the safety region on a particular side, and if so, don't allow entity transfers from that side
                for ett = 1 : length(Cell(i).Entities)
                    etx = [Cell(i).Entities(ett).x(1) Cell(i).Entities(ett).x(2) 0]; %set up for cross product
                    as = [Cell(i).side(sd).vertices(1,:) 0]; %set up for cross product
                    bs = [Cell(i).side(sd).vertices(2,:) 0]; %set up for cross product
                    aras = 0.5 * norm(cross(etx - as, etx - bs));  %get the area
                    dab = norm(Cell(i).side(sd).vertices(1,:) - Cell(i).side(sd).vertices(2,:), 2);
                    hts = aras * 2 / dab;
                    if (hts <= safeDistance)
                        Cell(i).side(sd).safe = 0; %if within the range, then it's not safe
                        break; % short circuit: exit for loop after detecting at least 1 entity which is unsafe
                    end
                    % todo: try to get method using linear inequalities for each side working, but this one works now
                end
            end
        end
    end
    
    for i = NF
        % SIGNALLING: i gives permission to j to move
        for j = Cell(i).nfnbrs
            % find the side corresponding to neighbor j (otherwise we don't have to do anything, since a side without a neighbor never has entities coming from it, nor going to it)
            sd = [];
            for stmp = 1 : Cell(i).num_sides
                if j == Cell(i).side(stmp).nbr
                    sd = stmp;
                    break;
                end
            end
            
            if isempty(sd)
                continue;
            end
        
            % let another side try to move (fairness)
            %if (~isempty(Cell(i).lastSideSignal) && Cell(i).lastSideSignal == j)
               %continue; % todo: may cause deadlocks? 2011-11-04
            %end

            % j has entities to move toward i of any color, the side joining i and j is safe
            % * Uses CellOld and Cell for j's Entities (to slow it down)
            % * Uses Cell for safety, since it was just recomputed
            % * Ensures entity types do not mix (i and j must have the same entity type variable)
            % * If cell i needs a lock to move, it has it
            c = Cell(i).etype;
            cj = Cell(j).etype;
            if ~isempty(Cell(i).side(sd).nbr) && Cell(i).side(sd).safe && Cell(i).side(sd).nbr == j && ~isempty(Cell(j).Entities)
                % color matching (could be anded with previous conditional, but for readability)
                % if cell i has entities, then they are the same color as cell j, whose next pointer of their shared color points at cell i
                % otherwise, cell i has no entities, so it can accept entities of any color from j
                % also: % i isn't resetting lock
                if ((~isempty(c) && c == cj && Cell(j).next( cj ) == i) || (isempty( c ) && Cell(j).next( cj ) == i))
                    if ~(Cell(i).detectNextFailed > 0 && k < Cell(i).detectNextFailed + opt_firstFullRound)
                    %if (isempty(Cell(j).path(cj).nlock) || (~isempty(Cell(j).path(cj).nlock) && Cell(j).path(cj).lock)) % locking
                    if (isempty(Cell(i).path(cj).nlock) || (~isempty(Cell(i).path(cj).nlock) && Cell(i).path(cj).lock)) % locking
                        if (~Cell(i).path(cj).rlock) % locking fairness; does use cj on purpose
                        % shouldn't need this in prev: || ~isempty(Cell(j).Entities)
                            Cell(i).side(sd).signal = j;
                            Cell(i).signal = j;
                            Cell(i).lastSideSignal = j;
                            if mode_debug
                                ['Cell ' num2str(i) ' signals to Cell ' num2str(j) ' on side ' num2str(sd)] % ' and Cell ' num2str(j) ' has next to Cell ' num2str(Cell(j).next( Cell(i).etype )) ' on side ' num2str(sd)]
                            end
                        end
                    end
                    end
                        %break;
                    %end
                end
            else
                Cell(i).lastSideSignal = [];
            end
        end
    end
     
    for i = NF
        % MOVING
        if ~isempty(Cell(i).Entities)
            tt = Cell(i).Entities( 1 ).color; % one color per cell, just grab the first entity's color
        else
            continue;
        end
        
        tp = Cell(i).next(tt); %the index of cell(i)'s next
        
        % todo: weird bug... should be prevented by signalling logic
        if tp == i
            continue;
        end
        
        % get side of next pointer; unused currently, was used when more than one color per cell allowed
        %for sd = 1: Cell(i).num_sides
        %    if Cell( tp ).side(sd).nbr == i
        %        sdj = sd;
        %        break;
        %    end
        %end
        
        % if the next cell is signalling to me, I can move toward it
        if (Cell( tp ).signal == i)
            %next line was used in previous, but appears unnecessary, could be used if more than one color allowed per cell: && Cell( tp ).side(sdj).signal == i)
            % find the side (face) in the direction of next
            for s = 1 : Cell(i).num_sides
                if Cell(i).side(s).nbr == tp
                    sd = s;
                end
            end

            if Cell(i).num_sides == 3
                % todo: move next two lines to triangulation creation
                vertCommon = intersect(roundn(Cell(tp).vertices,opt_decimal_tolerance), roundn(Cell(i).vertices,opt_decimal_tolerance), 'rows'); % common vertices between cells i and tp
                vertFarthest = setdiff(roundn(Cell(i).vertices,opt_decimal_tolerance), vertCommon, 'rows'); %the farthest vertices from the common edge
                Cell(i).ud = Cell(i).centroid - vertFarthest; % vector from the farthest vertice of Cell i through the centroid
            elseif Cell(i).num_sides == 4
                Cell(i).ud = Cell(i).side(sd).vectorMovement;
            else
                Cell(i).ud = Cell(i).side(sd).vectorMovement;
            end
            Cell(i).ud = Cell(i).ud / norm(Cell(i).ud, 2); % normalize the vector
            uds = (Cell(i).speed) * Cell(i).ud; % the position difference that the enitities need to move
            sj = 1;
            while sj <= length(Cell(i).Entities)
                % don't move entities twice (since we're using loops and not a perfectly synchronous model)
                if Cell(i).Entities(sj).moved
                    sj = sj + 1;
                    continue;
                end
                
                orig = Cell(i).Entities(sj).x; % original position
                temp = Cell(i).Entities(sj).x + uds; % new position
                if isempty(find(Cell(i).illegalRegion.A * temp' <= (Cell(i).illegalRegion.b + epsIllegal) == 0))
                    Cell(i).Entities(sj).x = temp;
                    Cell(i).Entities(sj).moved = 1;
                elseif tp == targets(tt) && tt == Cell(i).Entities(sj).color  %if the next cell is the target of the entity type's color, then remove the entity of that color (e.g., remove red entities on the red target, but not on the green target)
                    %pause
                    Cell(i).Entities = [Cell(i).Entities(1:sj - 1), Cell(i).Entities(sj + 1 : length(Cell(i).Entities))]; % remove the entity if the destination is target                    
                    sj = sj - 1; % gets added later for identity (but the max length will have decreased)
                    %Cell(tp).throughput(tt) = Cell(tp).throughput(tt) + 1.0; % counting the throughput
                    throughput(tt) = throughput(tt) + 1;
                    if min(size(Cell(i).Entities)) == 0 || isempty(Cell(i).Entities)
                        Cell(i).Entities = [];
                        Cell(i).lastEtype = Cell(i).etype; % last entity type on this cell will not be allowed in again until any other color requesting in gets in
                    end
                    
                    %if Cell(tp).path( Cell(i).etype ).lock %&& ((~isempty(intersect(sources, i)) && isempty(intersect(Cell(i).path( Cell(i).etype).int,i)))  || min(size(Cell(i).Entities)) == 0 || isempty(Cell(i).Entities))
                        %Cell(i).Entities
                        %Cell(tp).path( Cell(i).etype ).rlock = 1; % reset lock
                    %end
                else
                    %pause
                    od = 0;
                    st = length(Cell(tp).Entities); %get the original number of next cell's entities
                    Cell(tp).Entities(st + 1).x = temp; % add the transfered entity to the next cell
                    Cell(tp).Entities(st + 1).id = Cell(i).Entities(sj).id; % copy the index
                    Cell(tp).Entities(st + 1).color = Cell(i).Entities(sj).color; %add the transfered entity to the next cell
                    Cell(tp).Entities(st + 1).moved = 1;
                    Cell(i).Entities = [Cell(i).Entities(1:sj - 1), Cell(i).Entities((sj + 1): length(Cell(i).Entities))]; % remove the entity
                    sj = sj - 1; % DO NOT REFER TO sj AFTER THIS LINE
                    if min(size(Cell(i).Entities)) == 0 || isempty(Cell(i).Entities)
                        Cell(i).Entities = [];
                        Cell(i).lastEtype = Cell(i).etype; % last entity to be on this cell will not be allowed in again until any other color requesting in gets in
                    end
                    
                    %if Cell(tp).path( Cell(i).etype ).lock %&& ((~isempty(intersect(sources, i)) && isempty(intersect(Cell(i).path( Cell(i).etype).int,i)))  || min(size(Cell(i).Entities)) == 0 || isempty(Cell(i).Entities))
                        %Cell(i).Entities
                        %Cell(tp).path( Cell(i).etype ).rlock = 1; % reset lock
                    %end
                    
                    %'entity transfer'
                    a = orig; % a and b are two points on the movement line
                    b = temp;
                    [A,bv] = vert2con([a; b]);
                    
                    for nextSide = 1 : Cell(tp).num_sides
                        % have to sort, otherwise might not say equal due to different permutation
                        if isempty(find(sortrows(roundn(Cell(tp).side(nextSide).vertices, opt_decimal_tolerance)) == sortrows(roundn(Cell(i).side(sd).vertices, opt_decimal_tolerance)) == 0))
                            a = Cell(tp).illegalRegion.side(nextSide).vertices(1,:);
                            b = Cell(tp).illegalRegion.side(nextSide).vertices(2,:);
                        end
                    end

                    [As,bvs] = vert2con([a; b]);
                    Aceq = [A; As];
                    bceq = [bv; bvs];
                    
                    Cell(tp).Entities(st + 1).x = linprog(zeros(m,1),[],[],Aceq,bceq,[],[],[],opt_linprog_options)'; % point on the cell edge/boundary
                    vect = Cell(i).ud;
                end

                sj = sj + 1; % except here for increment
            end
        end
    end
    
    % sanity checks
    % 1. safety invariant check between all entities on each Cell
    % 2. all entities must be in the set of entities for the geographic
    %    position of the cell they reside on
    %    That is: an entity must be on the triangle it should be on
    % 3. check all entity ids to ensure no duplication
    for i = N
        fail_safety = 0;
        p = 1;
        while p <= length(Cell(i).Entities)
            % remove entities on targets of their color, if any managed to get there and not already be removed...
            if ~isempty(intersect(i, targets)) && Cell(i).color == Cell(i).Entities(p).color
                Cell(i).Entities = [Cell(i).Entities(1:p - 1) Cell(i).Entities(p + 1 : length(Cell(i).Entities))]; % remove the bad entity
                % don't change p, increment and decrement simultaneously
                continue;
            end
            
            % correct Entity set check
            % note use of epsIllegal: this is to avoid numerical rounding problems
            if ~isempty(find(Cell(i).illegalRegion.A * Cell(i).Entities(p).x' <= (Cell(i).illegalRegion.b + epsIllegal) == 0))
                ['PROPERTY VIOLATION: Entity is not on the correct Cell: it should be on Cell ', num2str(i), '; entity has id ', num2str(Cell(i).Entities(p).id), ' (x=', num2str(Cell(i).Entities(p).x(1)), ', y=', num2str(Cell(i).Entities(p).x(2)), ',local ids p=', num2str(p), ')']
                switch opt_badCell
                    case 0
                        % option 0: do nothing
                    case 1
                        % option 1: remove the bad entity
                        Cell(i).Entities = [Cell(i).Entities(1:p - 1) Cell(i).Entities(p + 1 : length(Cell(i).Entities))]; % remove the bad entity
                        p = p - 1;
                        if opt_display
                            call_plotSystem;
                        end
                
                    case 2
                        % option 2: reposition the bad entity
                        safe = 0;
                        while ~safe
                            % generate a random point in the source (convex combination of the vertices => always inside the convex set)
                            ccSides = rand(Cell(i).num_sides,1); % need a point from each side
                            while sum(ccSides(1: Cell(i).num_sides - 1)) >= 1
                                ccSides = rand(Cell(i).num_sides,1); % need a point from each side
                            end
                            ccSides( Cell(i).num_sides ) = 1 - sum(ccSides(1: Cell(i).num_sides - 1)); % ensures cc1 + cc2 + cc3 + ... = 1 and cc1,cc2,cc3,... > 0 (convex combination)
                            x = [0 0];
                            for side = 1 : Cell(i).num_sides
                                x = x + ccSides(side,:) * Cell(s).safeRegion.vertices(side,:);
                            end
                            for q = 1 : length(Cell(i).Entities)
                                if norm( Cell(i).Entities(p).x - x, 2) >= 2*d
                                    safe = 1;
                                else
                                    safe = 0;
                                    break;
                                end
                            end
                        end

                        % only reposition the entity if it's safe
                        if safe
                            Cell(i).Entities(p).x = x;
                        end
                end
                %pause;
            end
            
            for q = 1 : length(Cell(i).Entities)
                % duplicate entity check
                if p ~= q && Cell(i).Entities(p).id == Cell(i).Entities(q).id
                    ['ERROR: Entity duplication check failure on Cell ' num2str(i) ' between Entities with id ' num2str(Cell(i).Entities(p).id) ' (local ids p=' num2str(p) ', q=' num2str(q) ')']
                    if opt_display
                        call_plotSystem; % make sure we update the plot, otherwise it is inaccurate (previous round)
                    end
                    %pause
                end
                
                % safety invariant check (only display once)
                if p ~= q && ~fail_safety && norm(Cell(i).Entities(p).x - Cell(i).Entities(q).x, 2) < (2*L + rs)
                    ['PROPERTY VIOLATION: Safety invariant failure on Cell ' num2str(i)  ' between Entities with ids ' num2str(Cell(i).Entities(p).id) ' and ' num2str(Cell(i).Entities(q).id) ' (local ids p=' num2str(p) ', q=' num2str(q) ')']
                    fail_safety = 1;
                    if opt_display
                        call_plotSystem; % make sure we update the plot, otherwise it is inaccurate (previous round)
                    end
                    %pause
                end
            end
            
            %todo sanity checks
            % safety: see if there are two colors on any cell
            %
            % deadlock: see if there are two colors on any lockset
            
            p = p + 1;
        end
    end
    
    % reset moved flags
    for i = NF
        for p = 1 : length(Cell(i).Entities)
            Cell(i).Entities(p).moved = 0;
        end
    end
    
    % reset signals
    for i = NF
        for j = 1 : 3
            Cell(i).side(j).signal = i;
        end
    end
    
    CellOld = Cell; % copy for next round
    
    % PLOTTING
    if opt_display && mod(k,opt_display_update)==0
        clf; % clear the figure
        call_plotSystem; % call visualization routine
        if opt_waitplot
            pause(waitplot); % pause
        end

        % record each frame into the Movie object; indexing starts at k-1
        if opt_video
            %Movie(k-1) = getframe(gcf);
            frame = getframe;
            writeVideo(writerObj,frame);
        end
    end
    
    if mod(k,opt_display_progress) == 0
        ['System progressing, at round: ', num2str(k)]
        call_plotSystem;
        pause(waitplot);
    end
    
    %pause % manual update, press a key between each round
end

% save the video to file
if opt_video
    %movie2avi(Movie, 'factory.avi', 'compression', 'none', 'fps', 5);
    close(writerObj);
end

% for ti = 1 : NT
%      %Cell(targets(ti)).throughput = Cell(targets(ti)).throughput / k; % / MAX_ROUNDS
%      %target_tp(ti) = Cell(targets(ti)).throughput;
% end
throughput = throughput ./ k;

        simulations(iSim).run(iRun).throughput = throughput;
        simulations(iSim).run(iRun).entitiesPlaced = indexEntity;


    end

end