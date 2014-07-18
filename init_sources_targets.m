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