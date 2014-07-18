% associate neighbors with sides
% needs to be a separate loop: relies on Cell(j).vertices being constructed for j \in nbrs(i)
for i = 1 : N
    for s = 1 : Cell(i).num_sides
        Cell(i).side(s).safe = 1;   % side s starts safe (we run safety check before moving)
        
        for j = Cell(i).nbrs
            % delete bad neighbor indices
            if j > N || j <= 0
                Cell(i).nbrs = Cell(i).nbrs(Cell(i).nbrs ~= j); % remove neighbors with too large or too small indices
                continue;
            elseif isnan(j)
                Cell(i).nbrs = Cell(i).nbrs(~isnan(Cell(i).nbrs)); % remove NaN neighbors
                continue;
            end
            
            if size(intersect(roundn(Cell(i).side(s).vertices, opt_decimal_tolerance), roundn(Cell(j).vertices, opt_decimal_tolerance), 'rows')) == [2 2] % share two vertices => share an edge => adjacent triangles and this side is not on the boundary of the triangulation
                Cell(i).side(s).nbr = j;    % index of neighbor which touches cell i on side s
                Cell(i).side(s).boundary = 0;
            end
        end
    end
end