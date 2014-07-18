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