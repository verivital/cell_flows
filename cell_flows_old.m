
numpts = 25; % number of points to triangulate
L = 0.025; % entity radius
m = 2; % number of dimension (use this variable when necessary to refer to the dimensions)

pts = rand(numpts,m);
%pts = ndgrid(1:numpts,1:m) + rand(numpts,m)*0.01;
%pts = [0 0; 0 1; 1 0; 1 1; 2 0; 0 2; 2 2; 1 2; 2 1];
%pts = [0 0; 0 1; 1 0; 1 1; 2 0; 0 2; 2 2];

% todo: add constraints describing polygon being triangularized
dt = DelaunayTri(pts(:,1),pts(:,2));
inside = dt.inOutStatus();

% Construct a TriRep to represent the domain triangles.
tr = TriRep(dt(inside, :), dt.X);

% Construct a set of edges that join the circumcenters of neighboring
% triangles; the additional logic constructs a unique set of such edges.
N = size(tr,1); % number of cells/triangles
T = (1:N)';
neigh = tr.neighbors(); % set of all neighbors
cc = tr.circumcenters(); % set of all circumcenters
ic = tr.incenters(); % set of all incenters
idx1 = T < neigh(:,1); % modify indexing to make appropriate neighbor graph
idx2 = T < neigh(:,2);
idx3 = T < neigh(:,3);
neigh = [T(idx1) neigh(idx1,1); T(idx2) neigh(idx2,2); T(idx3) neigh(idx3,3)]';

targets = [1; 4];   % identifiers of the sinks (target cells)
sources = [9; 10];  % identifiers of the source cells: note, must be same length as targets, as each source will produce entities of the same type of the corresponding entry of the targets vector

MAX_ROUNDS = N + 1; % needs to be at least N (for routing to have necessarily stabilized), change as necessary

% Iterate over all N cells, initializing an array of structs of length N
%   * generate the neighbor graph for the cells using data structures generated
%     by the DelaunayTri command
%   * Note: cell is a keyword in Matlab, so make sure you use the word Cell to
%     when using this variable
for i = 1 : N
    Cell(i).nbrs = tr.neighbors(i);                 % copy neighbors from triangulation data
    Cell(i).dist = inf * ones(length(targets),1);   % infinite distance to every target initially
    Cell(i).next = i * ones(length(targets),1);     % point at oneself initially
    Cell(i).signal = [];                            % signal used by Cell i to indicate to a neighbor that it is safe to move in the direction of Cell i (i.e., no entities on Cell i along their common border)
    Cell(i).token = [];                             % token used to fairly signal each neighbor
    Cell(i).failed = 0;                             % initially cells are not faulty
    Cell(i).Entities = [];                          % set of entities on Cell i
    Cell(i).centroid = ic(i,:);                     % point in the plane of the center of the cell (for plotting, etc); this is a vector
    for tt = 1 : length(targets)
        Cell(i).epsPlot(tt,:) = [0.01*rand(1,1), 0.01*rand(1,1)];              % small spacing off the center for easier visualizing (see plotSystem)
    end
end

% Initialize all source cells
for i = 1 : length(sources)
    Cell(sources(i)).color = i;                             % set the "color" / type to be equal to the entry in the targets vector, i.e., 1, 2, 3, etc.
    Cell(sources(i)).Entities(1).x = Cell(sources(i)).centroid;      % start all the sources with a single entity (just for reference, eventually will need to place them)
    Cell(sources(i)).Entities(1).color = i;                 % set the color of the entity to be the color of the source
end

% Initialize all target cells to be 0 distance away from the target
for i = 1 : length(targets)
    Cell(targets(i)).dist(i) = 0;   % set the distance to 0 away from this color
    Cell(targets(i)).color = i;     % set the "color" / type to be equal to the entry in the targets vector, i.e., 1, 2, 3, etc.
end

waitplot = 0.05; % twentieth of a second between plots

CellOld = Cell; % copy the current cell information (models communication being delayed a round)
                % To model communications, from cell i, any access to
                % information of neighbors must be from the CellOld
                % variable, which represents what was communicated to i
                % from a neighbor at the last round.  For the current
                % round, update the Cell variable.

% main simulation loop, k is the round index
for k = 2 : MAX_ROUNDS
    % iterate over all the targets (create a different route for different color types; note that we are assuming there is a unique "color" or type corresponding to each target)
    for tt = 1 : length(targets)
        % iterate over all the cells
        for i = 1 : N
            % ROUTING
            % iterate over Cell i's neighbors
            for j = CellOld(i).nbrs
                if isnan(j) || j > N || j <= 0
                    continue;
                end
                
                % don't change sinks' pointers or distances
                if targets(tt) == i
                    continue;
                end

                if (CellOld(j).dist(tt) < CellOld(i).dist(tt))
                    Cell(i).dist(tt) = CellOld(j).dist(tt) + 1;
                    Cell(i).next(tt) = j;
                end
            end
            
            % SIGNALLING
            % todo, roughly going to be the same as the old version
            
            % MOVING
            if ~isempty(Cell(i).Entities)
                % roughly: for each entity of color tt, move in direction of next
            end
        end
    end
    
    CellOld = Cell; % copy for next round
    
    % PLOTTING
    clf; % clear the figure
    plotSystem; % call visualization routine
    pause(waitplot); % pause
end
