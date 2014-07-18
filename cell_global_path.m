
    
    
% initialize path variables
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