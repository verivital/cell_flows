% visualization routine
function out = plotSystem(Cell, targets, sources, failed, pts, N, L, rs, d)
    % drawing options
    opt_drawSafety = 0; % show the safety circles and regions
    opt_drawSafetySide = 0; % show if a side is safe or not
    opt_drawTransfer = 0; % show the transfer (aka illegal) circles and regions
    opt_drawIncircle = 0; % show the in-circle
    opt_drawNormalVector = 0; % show normal vectors to each face
    opt_drawMoveVector = 0; % show the movement vector (direction)
    opt_drawCellId = 1; % cell identifiers
    opt_drawEndId = 1; % source and target text
    opt_drawCentroid = 0; % show the center (centroid) of the cell
    opt_drawIncenter = 0; % show the incenter of the cell
    opt_drawNext = 0; % show the next pointer graph
    opt_drawPath = 1; % draw each color's path (not the same as route, from every source to target and from every cell with entities to target)
    opt_drawPathArrows = 1; % draw arrows on color's path
    opt_drawSourceTarget = 1; % draw the source and target filled in
    opt_entityId = 0; % show entity ids
    opt_drawLocks = 0; % show locks needed and obtained
    
    opt_drawFailed = 1; % 1 = show failed cells, 0 = don't
    opt_drawFailedText = 1; % 1 = show fail ids, 0 = don't
    
    axis off; % disable axes (paper images)

    faceColor = 'k';
    faceLineWidth = 3;
    pathLineWidth = 2;

    idFontSize = 12;
    endFontSize = 12;

    colors(1,:) = [0 1 1]; % cyan
    %colors(2,:) = [1 0 0]; % red
    %colors(3,:) = [0 1 0]; % green
    %colors(4,:) = [1 0 1]; % magenta
    colors(2,:) = [70 130 180]; % steel blue
    colors(3,:) = [205 92 92]; % indian red
    colors(4,:) = [102 205 170]; % med. aquamarine
    colors(5,:) = [221 160 221]; % plum
    colors(6,:) = [218 165 32]; % goldenrod
    colors = colors ./ 255;
    basecolors = size(colors,1);
    for tt = 1 : length(targets)
        r = 0.5*mod(tt,2);
        g = 0.33*mod(tt,3);
        b = 0.25*mod(tt,4);
        colors( basecolors + tt,:) = [r g b];
    end
    
    clf; % clear figure
    
    tree = [];
    for tt = 1 : length(targets)
        hold on;

        if opt_drawNext
            for i = 1 : N
                xn = [Cell(i).centroid(1); Cell(Cell(i).next(tt)).centroid(1)] + Cell(i).epsPlot(tt,1);  % make points representing a line from the center of a cell to its next pointer
                yn = [Cell(i).centroid(2); Cell(Cell(i).next(tt)).centroid(2)] + Cell(i).epsPlot(tt,2);
                if opt_drawPathArrows
                    a = [xn(1); yn(1)];
                    b = [xn(2); yn(2)];
                    u = (b - a);% / norm(b - a, 2);
                    quiver(xn(1), yn(1), u(1), u(2), 'Color', colors(mod(tt, length(colors)) + 1,:), 'LineWidth', pathLineWidth, 'MaxHeadSize', 2);           % plot the routing graph on the tesselation
                else
                    plot(xn, yn, 'Color', colors(mod(tt, length(colors)) + 1,:));           % plot the routing graph on the tesselation
                end

                % generate the routing tree for plotting, see plotRoutingTree and help treeplot
                if Cell(i).dist(tt) == 0
                    tree(i,tt) = 0;
                else
                    tree(i,tt) = Cell(i).next(tt);
                end
            end
        end

        % TODO: remove all this specific vertex order stuff, see poly2fv, poly2cw, poly2bool, etc.
        if opt_drawSourceTarget
            % fill in the target cells so they are clearly visible
            idxs = targets(tt);
            if Cell(idxs).num_sides == 3
                patch(Cell(idxs).vertices(:,1),Cell(idxs).vertices(:,2),colors(mod(tt, length(colors)) + 1,:));
            elseif Cell(idxs).num_sides == 4
                patch([Cell(idxs).vertices(1,1); Cell(idxs).vertices(3,1); Cell(idxs).vertices(4,1); Cell(idxs).vertices(2,1); Cell(idxs).vertices(1,1)], [Cell(idxs).vertices(1,2); Cell(idxs).vertices(3,2); Cell(idxs).vertices(4,2); Cell(idxs).vertices(2,2); Cell(idxs).vertices(1,2)], colors(mod(tt, length(colors)) + 1,:));
            end

            % fill in the source cells so they are clearly visible
            idxs = sources(tt);
            if Cell(idxs).num_sides == 3
                patch(Cell(idxs).vertices(:,1),Cell(idxs).vertices(:,2),colors(mod(tt, length(colors)) + 1,:));
            elseif Cell(idxs).num_sides == 4
                patch([Cell(idxs).vertices(1,1); Cell(idxs).vertices(3,1); Cell(idxs).vertices(4,1); Cell(idxs).vertices(2,1); Cell(idxs).vertices(1,1)], [Cell(idxs).vertices(1,2); Cell(idxs).vertices(3,2); Cell(idxs).vertices(4,2); Cell(idxs).vertices(2,2); Cell(idxs).vertices(1,2)], colors(mod(tt, length(colors)) + 1,:));
            end
        end
    end
    
    % fill in the failed cells
    if opt_drawFailed
        for i = 1 : length(failed)
            ix = failed(i);
            if ix <= N
                if Cell(ix).num_sides == 4
                    patch([Cell(ix).vertices(1,1); Cell(ix).vertices(3,1); Cell(ix).vertices(4,1); Cell(ix).vertices(2,1); Cell(ix).vertices(1,1)], [Cell(ix).vertices(1,2); Cell(ix).vertices(3,2); Cell(ix).vertices(4,2); Cell(ix).vertices(2,2); Cell(ix).vertices(1,2)], 'k');
                else
                    patch(Cell(ix).vertices(:,1), Cell(ix).vertices(:,2), 'k');
                end
            end
        end
    end

    % draw cell faces, show velocity vectors, etc.
    for i = 1 : N
        % show velocity vector if non-empty (bug checking)
        if opt_drawMoveVector && ~isempty(Cell(i).ud)
            quiver(Cell(i).centroid(1), Cell(i).centroid(2), Cell(i).ud(1), Cell(i).ud(2), 'k', 'MaxHeadSize', 1);
        end
        
        % only show non-faulty cells, or all cells if option enabled
        if ~Cell(i).failed || opt_drawFailed
            % draw the cell's faces / edges
            if Cell(i).num_sides == 3
                plot([Cell(i).vertices(:,1); Cell(i).vertices(1,1)], [Cell(i).vertices(:,2); Cell(i).vertices(1,2)], 'LineWidth', faceLineWidth, 'Color', faceColor);
            elseif Cell(i).num_sides == 4
                plot([Cell(i).vertices(1,1); Cell(i).vertices(3,1); Cell(i).vertices(4,1); Cell(i).vertices(2,1); Cell(i).vertices(1,1)], [Cell(i).vertices(1,2); Cell(i).vertices(3,2); Cell(i).vertices(4,2); Cell(i).vertices(2,2); Cell(i).vertices(1,2)], 'LineWidth', faceLineWidth, 'Color', faceColor);
            else
                % copies points around
                plot([Cell(i).vertices(:,1); Cell(i).vertices(1,1)],[Cell(i).vertices(:,2); Cell(i).vertices(1,2)], 'LineWidth', faceLineWidth, 'Color', faceColor);
            end
        end
        
        % show the safety subset of the cell
        % copy 1st point to make triangle
        if opt_drawSafety && Cell(i).safeRegion.R > 0
            if Cell(i).num_sides == 3
                plot([Cell(i).safeRegion.vertices(:,1); Cell(i).safeRegion.vertices(1,1)], [Cell(i).safeRegion.vertices(:,2); Cell(i).safeRegion.vertices(1,2)], 'r');
            elseif Cell(i).num_sides == 4
                plot([Cell(i).safeRegion.vertices(1,1); Cell(i).safeRegion.vertices(3,1); Cell(i).safeRegion.vertices(4,1); Cell(i).safeRegion.vertices(2,1); Cell(i).safeRegion.vertices(1,1)], [Cell(i).safeRegion.vertices(1,2); Cell(i).safeRegion.vertices(3,2); Cell(i).safeRegion.vertices(4,2); Cell(i).safeRegion.vertices(2,2); Cell(i).safeRegion.vertices(1,2)], 'r');
            else
                plot([Cell(i).safeRegion.vertices(:,1); Cell(i).safeRegion.vertices(1,1)],[Cell(i).safeRegion.vertices(:,2); Cell(i).safeRegion.vertices(1,2)], 'r');
            end
        end
        
        % show locks
        if opt_drawLocks
            for tt = 1 : length(targets)
                if targets(tt) == i || sources(tt) == i
                    tc = 'k';
                else
                    tc = colors(mod(tt, length(colors)) + 1,:);
                end
                
                if Cell(i).path(tt).lock || Cell(i).path(tt).rlock
                    text(Cell(i).centroid(1) + 0.1, Cell(i).centroid(2) + 0.1, ['l', num2str( Cell(i).path(tt).lock ), 'r', num2str( Cell(i).path(tt).rlock )], 'Color', tc);
                end
            end
        end
        
        % show the transfer subset of the cell
        if opt_drawTransfer && Cell(i).illegalRegion.R > 0 % if R <= 0, means the triangle will be inverted (i.e., the transfer set is the whole cell and the safe set is the whole cell)
            if Cell(i).num_sides == 3
                plot([Cell(i).illegalRegion.vertices(:,1); Cell(i).illegalRegion.vertices(1,1)], [Cell(i).illegalRegion.vertices(:,2); Cell(i).illegalRegion.vertices(1,2)], 'c');
            elseif Cell(i).num_sides == 4
                plot([Cell(i).illegalRegion.vertices(1,1); Cell(i).illegalRegion.vertices(3,1); Cell(i).illegalRegion.vertices(4,1); Cell(i).illegalRegion.vertices(2,1); Cell(i).illegalRegion.vertices(1,1)], [Cell(i).illegalRegion.vertices(1,2); Cell(i).illegalRegion.vertices(3,2); Cell(i).illegalRegion.vertices(4,2); Cell(i).illegalRegion.vertices(2,2); Cell(i).illegalRegion.vertices(1,2)], 'c');
            else
                % copies points around
                plot([Cell(i).illegalRegion.vertices(:,1); Cell(i).illegalRegion.vertices(1,1)],[Cell(i).illegalRegion.vertices(:,2); Cell(i).illegalRegion.vertices(1,2)], 'c');
            end
        end

        % show the centroid of the cell
        if opt_drawCentroid
            scatter(Cell(i).centroid(1), Cell(i).centroid(2), '+', 'k');
        end
        
        % show the incenter of the cell
        if opt_drawIncenter
            scatter(Cell(i).incenter(1), Cell(i).incenter(2), '+', 'r');
        end
    end

    % show which sides violate safety
    if opt_drawSafetySide
        % do this afterward so it's on top
        for i = 1 : N
            % indicate if side is unsafe
            for sd = 1 : Cell(i).num_sides
                if ~Cell(i).side(sd).safe
                    plot(Cell(i).side(sd).vertices(:,1), Cell(i).side(sd).vertices(:,2), 'm--', 'LineWidth', 4);
                    text(Cell(i).side(sd).midpoint(1) + Cell(i).epsPlot(1,1)*10, Cell(i).side(sd).midpoint(2) + Cell(i).epsPlot(1,2)*10, num2str(i));
                elseif Cell(i).side(sd).boundary
                    plot(Cell(i).side(sd).vertices(:,1), Cell(i).side(sd).vertices(:,2), 'k-.', 'LineWidth', 3);
                end
            end
        end
    end

    for i = 1 : N
        % plot normal vectors to sides (bug checking)
        if opt_drawNormalVector
            for sd = 1 : Cell(i).num_sides
                quiver(Cell(i).side(sd).midpoint(1), Cell(i).side(sd).midpoint(2), Cell(i).side(sd).vectorNormal(1)/5, Cell(i).side(sd).vectorNormal(2)/5, 'k');
            end
        end        
        
        if opt_drawPath
            for tt = 1 : length(targets)
                if ~isempty(intersect(Cell(i).path(tt).ids, i))
                    xn = [Cell(i).centroid(1); Cell(Cell(i).next(tt)).centroid(1)] + Cell(i).epsPlot(tt,1);  % make points representing a line from the center of a cell to its next pointer
                    yn = [Cell(i).centroid(2); Cell(Cell(i).next(tt)).centroid(2)] + Cell(i).epsPlot(tt,2);
                    %xn = [Cell(i).centroid(1); Cell(Cell(i).next(tt)).centroid(1)];  % make points representing a line from the center of a cell to its next pointer
                    %yn = [Cell(i).centroid(2); Cell(Cell(i).next(tt)).centroid(2)];
                    if opt_drawPathArrows
                        a = [xn(1); yn(1)];
                        b = [xn(2); yn(2)];
                        u = (b - a);% / norm(b - a, 2);
                        quiver(xn(1), yn(1), u(1), u(2), 'Color', colors(mod(tt, length(colors)) + 1,:), 'LineWidth', pathLineWidth, 'MaxHeadSize', 2);           % plot the routing graph on the tesselation
                        %'LineStyle', ':'
                    else
                        plot(xn, yn, 'Color', colors(mod(tt, length(colors)) + 1,:), 'LineWidth', pathLineWidth);           % plot the routing graph on the tesselation
                    end
                end
            end
        end
        
        % show in-circle
        if opt_drawIncircle
            rectangle('Curvature', [1 1], 'Position', [Cell(i).incenter(1)-Cell(i).Rp, Cell(i).incenter(2)-Cell(i).Rp, 2*Cell(i).Rp, 2*Cell(i).Rp], 'EdgeColor', 'm');
        end
        % shown in-circle corresponding to smaller safety-check region
        if opt_drawSafety && opt_drawIncircle
            if Cell(i).safeRegion.R > 0
                rectangle('Curvature', [1 1], 'Position', [Cell(i).incenter(1)-Cell(i).safeRegion.R, Cell(i).incenter(2)-Cell(i).safeRegion.R, 2*Cell(i).safeRegion.R, 2*Cell(i).safeRegion.R], 'EdgeColor', 'r');
            end
        end
        if opt_drawTransfer && opt_drawIncircle
            if Cell(i).illegalRegion.R > 0
                rectangle('Curvature', [1 1], 'Position', [Cell(i).incenter(1)-Cell(i).illegalRegion.R, Cell(i).incenter(2)-Cell(i).illegalRegion.R, 2*Cell(i).illegalRegion.R, 2*Cell(i).illegalRegion.R], 'EdgeColor', 'c');
            end
        end
    end

    % draw entities on top of everything
    for i = 1 : N
        for p = 1 : length(Cell(i).Entities)
            % draw a circle centered at Entities(p).x of radius l
            % Position: x y w h, where x, y are bottom left corner, and w, h are lengths to top right corners
            %rectangle('Curvature', [1 1], 'Position', [Cell(i).Entities(p).x(1)-d, Cell(i).Entities(p).x(2)-d, 2*d, 2*d], 'FaceColor', 'k');
            %rectangle('Curvature', [1 1], 'Position', [Cell(i).Entities(p).x(1)-L, Cell(i).Entities(p).x(2)-L, 2*L, 2*L], 'FaceColor', 'Color', colors(mod(Cell(i).Entities(p).color, length(colors)) + 1,:));

            [x,y,z] = cylinder(d,200);
            x = x + Cell(i).Entities(p).x(1);
            y = y + Cell(i).Entities(p).x(2);
            fill(x(1,:),y(1,:), 'k')
            
            [x,y,z] = cylinder(L,200);
            x = x + Cell(i).Entities(p).x(1);
            y = y + Cell(i).Entities(p).x(2);
            fill(x(1,:),y(1,:), colors(mod(Cell(i).Entities(p).color, length(colors)) + 1,:))
            
            
            if opt_entityId
                text(Cell(i).Entities(p).x(1), Cell(i).Entities(p).x(2), num2str(Cell(i).Entities(p).id));
            end
        end
        
        if Cell(i).failed
            %tcolor = 'w'; % white has some problem with axes off, use light gray
            tcolor = [0.94 0.94 0.94];
        else
            tcolor = 'k';
        end
        
        if ~Cell(i).failed || (opt_drawFailed && opt_drawFailedText)
            if opt_drawCellId && opt_drawEndId
                if isempty(intersect(sources, i)) && isempty(intersect(targets, i)) && opt_drawEndId
                    text(Cell(i).centroid(1),Cell(i).centroid(2),num2str(i), 'Color', tcolor, 'FontSize', idFontSize);   % display identifier on the corresponding cell
                end

                if ~isempty(intersect(targets, i)) && opt_drawEndId
                    %text(Cell(i).centroid(1),Cell(i).centroid(2), [num2str(i), sprintf('\n'), ' (Target)']);   % display target near corresponding source cell
                    text(Cell(i).centroid(1),Cell(i).centroid(2), [num2str(i),'_T'], 'Color', tcolor, 'FontSize', idFontSize);   % display identifier on the corresponding cell
                    %text(Cell(i).centroid(1) - 0.25, Cell(i).centroid(2) - 0.175, ['(Target)']);   % display target near corresponding source cell
                end

                if ~isempty(intersect(sources, i)) && opt_drawEndId
                    %text(Cell(i).centroid(1),Cell(i).centroid(2), [num2str(i), sprintf('\n'), ' (Source)']);   % display source near corresponding source cell
                    %text(Cell(i).centroid(1) - 0.25, Cell(i).centroid(2) - 0.175, ['(Source)']);   % display target near corresponding source cell
                    text(Cell(i).centroid(1),Cell(i).centroid(2), [num2str(i),'_S'], 'Color', tcolor, 'FontSize', idFontSize);   % display identifier on the corresponding cell
                end
            elseif opt_drawEndId
                if ~isempty(intersect(targets, i)) && opt_drawEndId
                    %text(Cell(i).centroid(1),Cell(i).centroid(2), [num2str(i), sprintf('\n'), ' (Target)']);   % display target near corresponding source cell
                    %text(Cell(i).centroid(1), Cell(i).centroid(2), ['Target'], 'FontSize', endFontSize);   % display target near corresponding source cell
                    text(Cell(i).centroid(1)-Cell(i).epsPlot(1),Cell(i).centroid(2)-Cell(i).epsPlot(2), [num2str(i),'_S'], 'Color', tcolor, 'FontSize', idFontSize);   % display identifier on the corresponding cell
                end

                if ~isempty(intersect(sources, i)) && opt_drawEndId
                    %text(Cell(i).centroid(1),Cell(i).centroid(2), [num2str(i), sprintf('\n'), ' (Source)']);   % display source near corresponding source cell
                    %text(Cell(i).centroid(1), Cell(i).centroid(2), ['Source'], 'FontSize', endFontSize);   % display target near corresponding source cell
                    text(Cell(i).centroid(1)-Cell(i).epsPlot(1),Cell(i).centroid(2)-Cell(i).epsPlot(2), [num2str(i),'_S'], 'Color', tcolor, 'FontSize', idFontSize);   % display identifier on the corresponding cell
                end
            end
        end
    end

    axis equal;
    xlabel('x');
    ylabel('y');
    %pause
end