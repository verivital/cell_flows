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