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