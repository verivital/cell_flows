% compute the path
for i = NF
    for tt = 1 : NT
        if opt_routing_disjoint == 0
            % inductively defined as:
            % path_0: if i has entities or is a source
            % path_{i}: if i is in the ids of this color, add its' next

            % ensure routes have stabilized (sort of...)
            %if CellOld(i).next ~= Cell(i).next
            %    break; % don't do anyone if any have changed
            %end

            if ((~isempty(Cell(i).Entities) && Cell(i).etype == tt) || sources(tt) == i || ~isempty(intersect(Cell(i).path(tt).ids, i))) && Cell(i).next(tt) ~= i
                Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; i; Cell(i).next(tt)]);
            end

            % todo: after a long time from the last failure, we may be able to remove cells from
            % being on the path, e.g., after there are only entities on
            % cells between the source and path, we can remove the old ones

            % aggregate paths of all neighbors (gossip)
            for j = Cell(i).nbrs
                Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; CellOld(j).path(tt).ids]); % note use of old
            end
        else
            restart = 1;
            % reversal version
            while (restart)
                restart = 0;
                if ((~isempty(Cell(i).Entities) && Cell(i).etype == tt) || sources(tt) == i || ~isempty(intersect(Cell(i).path(tt).ids, i))) && Cell(i).next(tt) ~= i
                    Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; i; Cell(i).next(tt)]);
                end

                for tj = 1 : NT
                    if tt ~= tj
                        % path for tj subseteq tt
                        common = intersect(Cell(i).path(tt).ids,Cell(i).path(tj).ids);
                        if ~isempty(common) && length(common) == length(Cell(i).path(tj).ids) && isempty(find((common == Cell(i).path(tj).ids) == 0))
                            %'subset'
                            %for rev = 1 : length(Cell(i).path(tt).ids)
                            %    ri = Cell(i).path(tt).ids(rev);
                            %    tmp = ceil(length(Cell(ri).nfnbrs)*rand(1,1));
                            %    Cell(ri).next(tt) = Cell(ri).nfnbrs(tmp);
                            %end
                        end
                    end
                end

                % todo: after a long time from the last failure, we may be able to remove cells from
                % being on the path, e.g., after there are only entities on
                % cells between the source and path, we can remove the old ones

            end
            % aggregate paths of all neighbors (gossip)
            for j = Cell(i).nbrs
                Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; CellOld(j).path(tt).ids]); % note use of old
            end
        end

%             % a-star version
%             restart = 1;
%             retry = 0;
%             nextToTry = Cell(i).next(tt);
%             while restart
%                 %Cell(i).path(tt).ids = CellOld(i).path(tt).ids;
%                 restart = 0;
%                 if ((~isempty(Cell(i).Entities) && Cell(i).etype == tt) || sources(tt) == i || ~isempty(intersect(Cell(i).path(tt).ids, i))) && nextToTry ~= i && isempty(intersect(Cell(i).path(tt).nextTried, nextToTry))
%                     if retry
%                         Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; nextToTry]);
%                     else
%                         Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; i; nextToTry]);
%                     end
%                 end
% 
%                 % todo: after a long time from the last failure, we may be able to remove cells from
%                 % being on the path, e.g., after there are only entities on
%                 % cells between the source and path, we can remove the old ones
%                  
%                 Cell(i).path(tt).int = []; % reset to start empty
%                  
%                 for tj = 1 : NT
%                     if tt == tj
%                        continue
%                     end
%                     pint = unique(intersect( Cell(i).path(tt).ids, Cell(i).path(tj).ids ));
%                     Cell(i).path(tt).int = unique([Cell(i).path(tt).int; pint]);
%                 end
%  
%                 if ~isempty(Cell(i).path(tt).int) && ~isempty(intersect(Cell(i).path(tt).int, i))
%                     'retry'
%                     i
%                     tt
%                     
%                     Cell(i).path(tt).nextTried = unique([Cell(i).path(tt).nextTried; nextToTry]);
%                     
%                     tmpnbrs = Cell(i).nfnbrs;
% %                     tmpnbrs = [];
% %                     for j = Cell(i).nfnbrs
% %                         %if Cell(j).next(tt) ~= i
% %                         if ~isempty(intersect(Cell(i).path(tt).ids, j))
% %                             tmpnbrs = [tmpnbrs; j];
% %                         end
% %                     end
% 
%                     tmp = setdiff(tmpnbrs, Cell(i).path(tt).nextTried);
%                     if isempty(tmp) % have to quit
%                         restart = 0;
%                         if i == 17 || i == 12 || i == 22 || i == 7 || i == 2
%                             Cell(i).path(tt).nextTried
%                             %pause
%                         end
%                     else
%                         len = length(tmp);
%                         Cell(i).path(tt).ids = setdiff(Cell(i).path(tt).ids, [nextToTry]);
%                         %Cell(i).path(tt).ids = [i];
%                         %nextToTry = tmp(1);
%                         nextToTry = tmp( ceil(len*rand(1,1)) );
%                         
%                         if i == 17 || i == 12 || i == 22 || i == 7 || i == 2
%                             'trying next'
%                             nextToTry
%                             %pause
%                         end
% 
%                         restart = 1;
%                         retry = 1;
%                     end
%                 end
%             end
%             
%             % aggregate paths of all neighbors (gossip)
%             for j = Cell(i).nfnbrs
%                 
%                 % if i isn't in it's own path, we need to not gossip it around
%                 %if isempty(intersect(Cell(i).path(tt).ids, i)) && ~isempty(intersect(CellOld(j).path(tt).ids, i))
%                 %    Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; CellOld(j).path(tt).ids]); % note use of old
%                 %    Cell(i).path(tt).ids = setdiff(Cell(i).path(tt).ids, i);
%                 %else
%                 %    Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; CellOld(j).path(tt).ids]); % note use of old
%                 %end
%                 Cell(i).path(tt).ids = unique([Cell(i).path(tt).ids; CellOld(j).path(tt).ids]); % note use of old
%             end
    end
end
