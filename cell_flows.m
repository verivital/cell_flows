%function target_tp = cell_flows(targets, sources)

% load all possible simulation scenarios (particular environments, sources,
% targets, etc. that have been saved, such as snub square tesselations, square 
% tesselations, parallelogram tesselations, etc.)
init_scenarios;

if ~opt_single
    min_sim = 1;
    max_sim = length(simulations);
else
    min_sim = opt_single;
    max_sim = opt_single;
    min_run = opt_run;
    max_run = opt_run;
end

% main simulation loop over different scenarios to simulate (for allowing
% batch processing and throughput statistics creations)
for iSim = min_sim : max_sim
    if ~opt_run
        min_run = 1;
        max_run = length(simulations(iSim).run);
    end
    
    % minor simulation loop over different subscenarios for a given
    % environment, e.g., different sources/targets for the same tesselation
    for iRun = min_run : max_run
        sources = simulations(iSim).run(iRun).sources;
        targets = simulations(iSim).run(iRun).targets;
        opt_tesselation_type = simulations(iSim).environment;
        try
            failed = simulations(iSim).run(iRun).failed;
        catch
            failed = [];
        end
        
        try % need to be separate blocks, each can fail independently
            faildynamic = simulations(iSim).run(iRun).faildynamic;
        catch
            faildynamic = [];
        end
        
        opt_triangulation = 0; % 0 = don't use delauny triangulation (use a uniformly generated set of vertices instead, for throughput simulations)
        
        placedColors = [];
        lockedColors = [];

        L = 0.1;
        %L = 0.05;       % entity radius; for throughput
        
        if opt_tesselation_type == OPT_TESSELATION_TRI || opt_tesselation_type == OPT_TESSELATION_EQTRI_LINE || opt_tesselation_type == OPT_TESSELATION_RTRI_LINE
            epsPlot = 0.01;
            L = L / 2;
        else
            epsPlot = 0.25;
        end
        epsIllegal = 0.0001;
        
        rs = L / 10;    % safety radius
        d = L + rs; 
        m = 2;          % number of dimension (use this variable when necessary to refer to the dimensions)
        F = 0;          % number of cells to fail
        safeDistance = 3 * L + rs;

        % non-batched options
        NT = 1; % number of targets
        NS = NT; % number of sources

        MAX_ROUNDS = simulations(iSim).maxRounds;

        % options (potentially per-scenario)
        mode_debug = 0;
        opt_randic = 0; % random initial condition?
        opt_fixedic = 1; % case number for fixed initial conditions (see switch/case statement below); only applies if opt_randic = 0
        opt_center_entity = 0; % 1 = place one entity at center of cell at a time, 0 = place (about) numToTryPlacing
        %numToTryPlacing = ceil(1/L)/2; % todo: generalize and find side length
        numToTryPlacing = 20;
        opt_roundsToPlace = 1; % try to place entities every opt_roundsToPlace rounds
        opt_display = 1;        % don't draw system evolution if 0, draw if 1
        %opt_display_update = 10;
        opt_display_update = 1;
        opt_display_progress = 100;  % number of rounds to display a text update on progress if we aren't drawing the system evolution
        opt_decimal_tolerance = -4; % decimal tolerance (10^opt_decimal_tolerance)
        opt_badCell = 1;
        opt_singleEntity = 0;
        opt_video = 1; % record video

        if opt_video
            writerObj = VideoWriter('factory.avi');
            writerObj.FrameRate = 5;
            open(writerObj);
        end

        opt_linprog_options=optimset('Display', 'off'); % used to supress messages from linear programming solver

        opt_routing_disjoint = 0; % choose routing mode: we may in the future want to implement the disjoint paths algorithm, but for now will not

        gcf;
        %set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        axis equal;
        set(gcf, 'Position', [1 1 800 800]); % Maximize figure.

        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');

        % set failed cells manually
        failround = 75;
        %failround = -1;
        opt_fail_random = 0;
        opt_faildynamic = 1; % option to allow failures occurring while running

        % initialize environment and partition thereof
        init_cell_partitions;

        % must be defined after diam (created during environment and partition initialization)
        opt_firstFullRound = 2*diam; % first full round to place entities, etc.; in general needs to be 2*N (worst-case diameter)

        % set up neighbor graph for triangulated environments
        if opt_triangulation
            T = (1:N)';
            neigh = tr.neighbors(); % set of all neighbors
            cc = tr.circumcenters(); % set of all circumcenters
            [ic,icr] = tr.incenters(); % set of all incenters with corresponding radii (inscribed circle center and radius)
            idx1 = T < neigh(:,1); % modify indexing to make appropriate neighbor graph
            idx2 = T < neigh(:,2);
            idx3 = T < neigh(:,3);
            neigh = [T(idx1) neigh(idx1,1); T(idx2) neigh(idx2,2); T(idx3) neigh(idx3,3)]';
        end

        % set sources and targets
        init_sources_targets;
        % initialize other cell variables
        init_cell_vars;
        % initialize path data structures
        init_paths;
        % initialize and set failure data structures
        init_failures;
        % error checking for numbers of sources, targets, failures, etc.
        check_sources_targets_failures;

        % initialize neighbors
        find_neighbors;

        % initialize side-to-neighbor data structures
        find_sides;

        % Initialize all source cells
        for i = 1 : length(sources)
            Cell(sources(i)).color = i;	% set the "color" / type to be equal to the entry in the targets vector, i.e., 1, 2, 3, etc.
        end

        % Initialize all target cells to be 0 distance away from the target
        for i = 1 : length(targets)
            Cell(targets(i)).dist(i) = 0;   % set the distance to 0 away from this color
            Cell(targets(i)).color = i;     % set the "color" / type to be equal to the entry in the targets vector, i.e., 1, 2, 3, etc.
        end

        opt_waitplot = 1; % on = 1, off = 0
        waitplot = 0.000001; % wait (in seconds) between plots

        indexEntity = 1; % global index for entities (to see if we duplicate, compare them between rounds, etc.)

        CellOld = Cell; % copy the current cell information (models communication being delayed a round)
                        % To model communications, from cell i, any access to
                        % information of neighbors must be from the CellOld
                        % variable, which represents what was communicated to i
                        % from a neighbor at the last round.  For the current
                        % round, update the Cell variable.
        locks = [];
        overlap = 0;
        startOver = ones(NT,1);

        % main simulation loop, k is the round index
        for k = 2 : MAX_ROUNDS
            NF = [1:N];
            % updated failed / non-failed cells
            for i = 1 : N
                if Cell(i).failed
                    NF = setdiff(NF, i);
                end

                Cell(i).nfnbrs = Cell(i).nbrs;
                % could use NF, but then would need another loop
                for j = Cell(i).nbrs
                    if Cell(j).failed
                        Cell(i).nfnbrs = setdiff(Cell(i).nfnbrs, j);
                    end
                end
            end

            for ti = 1 : NT
                if isempty(path(ti).int)
                    placedColors = setdiff(placedColors, ti);
                end
            end

            % compute previous pointers in the routing path
            cell_subroutine_prev;

            % fail cells dynamically
            cell_global_fail_dynamic;

% old version of locking
%     for ti = 1 : NT
%         for i = NF
%             if sources(ti) == i && ~isempty(intersect(placedColors, ti)) % this color already put stuff on the source, block it
%                 %Cell(i).path(ti).rlock = 1;
%             else
%                 %Cell(i).path(ti).rlock = isempty(Cell(i).Entities);
%             end
%         end
%     end
%     
%     for ti = 1 : NT
%         for i = NF
%             if ~(sources(ti) == i && ~isempty(intersect(placedColors, ti))) % this color already put stuff on the source, block it
%                 for j = Cell(i).nfnbrs
%                     %Cell(i).path(ti).rlock = Cell(i).path(ti).rlock & Cell(j).path(ti).rlock;
%                 end
%             end
%         end
%     end


            % compute global paths (emulating gossip over all cells)
            cell_global_path;

            % compute global locks over paths (emulating distributed mutual exclusion
            % algorithm over all cells)
            cell_global_locks;

            % try to place entities safely and without violating assumptions needed
            % for liveness
            % NOTE: must be performed after computing paths and locks, since this
            % is how we enforce the assumptions
            cell_place_entities;

            % color cells based on entity types
            for i = NF
                Cell(i).etype = [];
                for p = 1 : length(Cell(i).Entities)
                    if isempty(intersect(Cell(i).etype, Cell(i).Entities(p).color))
                        Cell(i).etype = [Cell(i).etype Cell(i).Entities(p).color];
                    end
                end
                Cell(i).etype = unique(Cell(i).etype); % ensure unique

                % if a cell doesn't have a previous entity type yet, go ahead and update
                % otherwise, it is only ever updated by the move function when the last entity leaves the cell
                if isempty(Cell(i).lastEtype)
                    Cell(i).lastEtype = Cell(i).etype;
                end

                for tt = 1 : NT
                    % see if next pointer is failed, detect failure
                    np = Cell(i).next(tt);
                    if isinf(Cell( np ).dist(tt)) && k >= opt_firstFullRound
                        Cell(i).detectNextFailed = k;
                    end
                end
            end

            % ROUTING: compute the routing graphs from all cells to targets
            cell_subroutine_route;

            % LOCKING: compute the paths from all cells with entities to targets
            cell_subroutine_path;

            % LOCKING: compute the locks for all paths (from all cells with entities to targets)
            cell_subroutine_locks;

            % SIGNALING: compute the signals for all cells and paths
            cell_subroutine_signal;

            % MOVING: move entities on cells it is safe to do so for (based on received signal)
            cell_subroutine_move;

            % runtime monitoring for safety and liveness properties
            cell_runtime_monitor;

            % reset moved flags for all entities
            for i = NF
                for p = 1 : length(Cell(i).Entities)
                    Cell(i).Entities(p).moved = 0;
                end
            end

            % reset signals for all cells and sides
            for i = NF
                for j = 1 : length(Cell(i).side)
                    Cell(i).side(j).signal = i;
                end
            end

            CellOld = Cell; % copy for next round

            % PLOTTING
            if opt_display && mod(k,opt_display_update)==0
                clf; % clear the figure
                call_plotSystem; % call visualization routine
                if opt_waitplot
                    pause(waitplot); % pause
                end

                % record each frame into the Movie object; indexing starts at k-1
                if opt_video
                    %Movie(k-1) = getframe(gcf);
                    frame = getframe;
                    writeVideo(writerObj,frame);
                end
            end

            if mod(k,opt_display_progress) == 0
                ['System progressing, at round: ', num2str(k)]
                call_plotSystem;
                pause(waitplot);
            end

            %pause % manual update, press a key between each round
        end % end main time-step / round loop

        % save the video to file
        if opt_video
            %movie2avi(Movie, 'factory.avi', 'compression', 'none', 'fps', 5);
            close(writerObj);
        end

        % for ti = 1 : NT
        %      %Cell(targets(ti)).throughput = Cell(targets(ti)).throughput / k; % / MAX_ROUNDS
        %      %target_tp(ti) = Cell(targets(ti)).throughput;
        % end
        throughput = throughput ./ k;

        simulations(iSim).run(iRun).throughput = throughput;
        simulations(iSim).run(iRun).entitiesPlaced = indexEntity;

    end
end