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