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