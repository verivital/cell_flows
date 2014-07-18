% find neighbors (inefficient all to all compare, could be done smarter)
for i = 1 : N
    for j = 1 : N
        % cells are neighbors iff they share an edge
        if i ~= j && size(intersect(roundn(Cell(i).vertices,opt_decimal_tolerance), roundn(Cell(j).vertices,opt_decimal_tolerance),'rows'),1) >= 2
            Cell(i).nbrs = unique([Cell(i).nbrs, j]);
            %Cell(j).nbrs = [Cell(j).nbrs, i];
        end
    end
    Cell(i).nfnbrs = Cell(i).nbrs;
end