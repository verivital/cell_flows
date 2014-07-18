% iterate over all the targets (create a different route for different color types; 
    % note that we are assuming there is a unique "color" or type corresponding to each target)
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