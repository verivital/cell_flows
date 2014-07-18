for i = NF
    % MOVING
    if ~isempty(Cell(i).Entities)
        tt = Cell(i).Entities( 1 ).color; % one color per cell, just grab the first entity's color
    else
        continue;
    end

    tp = Cell(i).next(tt); %the index of cell(i)'s next

    % todo: weird bug... should be prevented by signalling logic
    if tp == i
        continue;
    end

    % get side of next pointer; unused currently, was used when more than one color per cell allowed
    %for sd = 1: Cell(i).num_sides
    %    if Cell( tp ).side(sd).nbr == i
    %        sdj = sd;
    %        break;
    %    end
    %end

    % if the next cell is signalling to me, I can move toward it
    if (Cell( tp ).signal == i)
        %next line was used in previous, but appears unnecessary, could be used if more than one color allowed per cell: && Cell( tp ).side(sdj).signal == i)
        % find the side (face) in the direction of next
        for s = 1 : Cell(i).num_sides
            if Cell(i).side(s).nbr == tp
                sd = s;
            end
        end

        if Cell(i).num_sides == 3
            % todo: move next two lines to triangulation creation
            vertCommon = intersect(roundn(Cell(tp).vertices,opt_decimal_tolerance), roundn(Cell(i).vertices,opt_decimal_tolerance), 'rows'); % common vertices between cells i and tp
            vertFarthest = setdiff(roundn(Cell(i).vertices,opt_decimal_tolerance), vertCommon, 'rows'); %the farthest vertices from the common edge
            Cell(i).ud = Cell(i).centroid - vertFarthest; % vector from the farthest vertice of Cell i through the centroid
        elseif Cell(i).num_sides == 4
            Cell(i).ud = Cell(i).side(sd).vectorMovement;
        else
            Cell(i).ud = Cell(i).side(sd).vectorMovement;
        end
        Cell(i).ud = Cell(i).ud / norm(Cell(i).ud, 2); % normalize the vector
        uds = (Cell(i).speed) * Cell(i).ud; % the position difference that the enitities need to move
        sj = 1;
        while sj <= length(Cell(i).Entities)
            % don't move entities twice (since we're using loops and not a perfectly synchronous model)
            if Cell(i).Entities(sj).moved
                sj = sj + 1;
                continue;
            end

            orig = Cell(i).Entities(sj).x; % original position
            temp = Cell(i).Entities(sj).x + uds; % new position
            if isempty(find(Cell(i).illegalRegion.A * temp' <= (Cell(i).illegalRegion.b + epsIllegal) == 0))
                Cell(i).Entities(sj).x = temp;
                Cell(i).Entities(sj).moved = 1;
            elseif tp == targets(tt) && tt == Cell(i).Entities(sj).color  %if the next cell is the target of the entity type's color, then remove the entity of that color (e.g., remove red entities on the red target, but not on the green target)
                %pause
                Cell(i).Entities = [Cell(i).Entities(1:sj - 1), Cell(i).Entities(sj + 1 : length(Cell(i).Entities))]; % remove the entity if the destination is target                    
                sj = sj - 1; % gets added later for identity (but the max length will have decreased)
                %Cell(tp).throughput(tt) = Cell(tp).throughput(tt) + 1.0; % counting the throughput
                throughput(tt) = throughput(tt) + 1;
                if min(size(Cell(i).Entities)) == 0 || isempty(Cell(i).Entities)
                    Cell(i).Entities = [];
                    Cell(i).lastEtype = Cell(i).etype; % last entity type on this cell will not be allowed in again until any other color requesting in gets in
                end

                %if Cell(tp).path( Cell(i).etype ).lock %&& ((~isempty(intersect(sources, i)) && isempty(intersect(Cell(i).path( Cell(i).etype).int,i)))  || min(size(Cell(i).Entities)) == 0 || isempty(Cell(i).Entities))
                    %Cell(i).Entities
                    %Cell(tp).path( Cell(i).etype ).rlock = 1; % reset lock
                %end
            else
                %pause
                od = 0;
                st = length(Cell(tp).Entities); %get the original number of next cell's entities
                Cell(tp).Entities(st + 1).x = temp; % add the transfered entity to the next cell
                Cell(tp).Entities(st + 1).id = Cell(i).Entities(sj).id; % copy the index
                Cell(tp).Entities(st + 1).color = Cell(i).Entities(sj).color; %add the transfered entity to the next cell
                Cell(tp).Entities(st + 1).moved = 1;
                Cell(i).Entities = [Cell(i).Entities(1:sj - 1), Cell(i).Entities((sj + 1): length(Cell(i).Entities))]; % remove the entity
                sj = sj - 1; % DO NOT REFER TO sj AFTER THIS LINE
                if min(size(Cell(i).Entities)) == 0 || isempty(Cell(i).Entities)
                    Cell(i).Entities = [];
                    Cell(i).lastEtype = Cell(i).etype; % last entity to be on this cell will not be allowed in again until any other color requesting in gets in
                end

                %if Cell(tp).path( Cell(i).etype ).lock %&& ((~isempty(intersect(sources, i)) && isempty(intersect(Cell(i).path( Cell(i).etype).int,i)))  || min(size(Cell(i).Entities)) == 0 || isempty(Cell(i).Entities))
                    %Cell(i).Entities
                    %Cell(tp).path( Cell(i).etype ).rlock = 1; % reset lock
                %end

                %'entity transfer'
                a = orig; % a and b are two points on the movement line
                b = temp;
                [A,bv] = vert2con([a; b]);

                for nextSide = 1 : Cell(tp).num_sides
                    % have to sort, otherwise might not say equal due to different permutation
                    if isempty(find(sortrows(roundn(Cell(tp).side(nextSide).vertices, opt_decimal_tolerance)) == sortrows(roundn(Cell(i).side(sd).vertices, opt_decimal_tolerance)) == 0))
                        a = Cell(tp).illegalRegion.side(nextSide).vertices(1,:);
                        b = Cell(tp).illegalRegion.side(nextSide).vertices(2,:);
                    end
                end

                [As,bvs] = vert2con([a; b]);
                Aceq = [A; As];
                bceq = [bv; bvs];

                Cell(tp).Entities(st + 1).x = linprog(zeros(m,1),[],[],Aceq,bceq,[],[],[],opt_linprog_options)'; % point on the cell edge/boundary
                vect = Cell(i).ud;
            end

            sj = sj + 1; % except here for increment
        end
    end
end