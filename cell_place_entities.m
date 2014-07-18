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