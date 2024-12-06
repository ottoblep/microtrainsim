function [traj, events] = constructTrajectory(network, params, solution, initial_position, initial_speed, planned_stops)
    %% Constructs a single train trajectory on the graph
    % Also iteratively enforces scheduled edges, dead ends and speed limits by adjusting speed targets
    % Simulation is based on edge change events then trajectory is filled with continuous positions

    % Target velocity vars are a sparse representation of the target velocity curves.
    % They are made up of support points (time-value pairs) that are then interpolated during trajectory construction.
    % This reduces the variable count and produces less erratic movement.
    % Switch directions represent the choice of next edge when hitting a vertex 

    % solution dimensions ( 2 * n_v_target_vars + n_switch_vars )
    % solution values (v_target_time 0-1, v_target_value 0-1, switch_direction 0-1)
    % solution serialization (v_target_timesteps, v_target_values, sw_direction)

    % initial position dimensions (3)
    % initial position values (edge 0-n, position on edge 0-1, train orientation on edge -1, 1)

    % initial speed dimensions (1)
    % initial speed values (speed: -max_speed-max_speed)

    % planned stops dimenstions (n, 2)
    % planned stops values (edge 0-n, arrival_time_fraction_of_timesteps 0-1)

    % events (n_edge_changes, 2)
    % events (timestep, edge 1-n, position on edge 0-1, train orientation on edge -1, 1, speed -inf, inf)

    % trajectory dimensions (timestep, 3)
    % trajectory values (edge 1-n, position on edge 0-1, train orientation on edge -1, 1, speed -inf, inf)

    % Aquire switch directions from solution
    switch_directions = solution(2 * params.n_v_target_vars + 1:2 * params.n_v_target_vars + params.n_switch_vars);

    % Aquire speed target points from solution
    v_target_timesteps = ceil(solution(1:params.n_v_target_vars) * params.n_timesteps);
    [v_target_timesteps , unique_idxs, ~] = unique(v_target_timesteps, 'stable');
    v_target_values = solution(params.n_v_target_vars + 1:2 * params.n_v_target_vars);
    v_target_values = (2 * v_target_values(unique_idxs) - 1) * params.max_speed;
    v_targets = [v_target_timesteps' v_target_values'];

    % Events is the first train state for each new edge
    % Trajectory is the train state at each timestep

    % p(n) = p(n-1) + v(n-1)
    % v(n) = v(n-1) + max(verr(n-1), a_max)

    events(1,:) = [1, initial_position(1), initial_position(2), initial_position(3), initial_speed];
    traj(1,:) = events(1, 2:5); 

    abort = 0;
    while true
        assert(all(all(events(:, 1:2) ~= 0)));

        if (events(end, 5) > network.speed_limits(events(end,2)) || events(end, 5) > params.max_speed)
            error("Speed limit could not be satisfied.");
        end

        % Adjust copy of speed targets for current edge speed limit
        % Assumptions:
        %   - Braking can only shift edge traversals backwards
        %   - When jumping back multiple transitions these are inside the braking curve (fixed speed target)
        %   - Adding more braking in this way can never exceed a previous speed limit
        v_targets_working_set = v_targets;
        edge_target_idxs = find(v_targets(:,1) >= events(end,1));
        v_targets_working_set(edge_target_idxs, 2) = min(v_targets(edge_target_idxs, 2), network.speed_limits(events(end, 2)));

        % Find next edge exit
        [edge_transition edge_trajectory edge_speeds] = simulateEdge(network, params, events(end, :), v_targets_working_set);

        if isempty(edge_transition)
            % Finish simulation
            traj(events(end,1):params.n_timesteps, 1) = events(end, 2);
            traj(events(end,1):params.n_timesteps, 2) = edge_trajectory;
            traj(events(end,1):params.n_timesteps, 3) = events(end, 4);
            traj(events(end,1):params.n_timesteps, 4) = edge_speeds;
            break;
        end

        viable_next_edges = network.adjacent_edge_list{edge_transition.traversed_node};
        viable_next_edges = viable_next_edges(viable_next_edges~=events(end,1));

        % Further modify speed targets in case of constraint violations

        v_targets_modified = true;
        global_speeds = cat(1, traj(1:(events(end,1) - 1), 4), edge_speeds', edge_transition.speed);

        % Check for dead end
        if isempty(viable_next_edges)
            % Stop train until speed target is in the other direction
            reverse_speed_target_idxs = (v_targets_working_set(:,1) >= edge_transition.timestep & sign(v_targets_working_set(:, 2)) ~= sign(edge_speeds(end)));
            departure_timestep = min(v_targets_working_set(reverse_speed_target_idxs, 1));
            [v_targets_working_set start_braking_timestep] = addBraking(params, global_speeds, v_targets_working_set, edge_transition, false, departure_timestep, 0);
        else
            % Decide next edge
            next_edge_selection = 1 + round(switch_directions(size(events, 1)) * (length(viable_next_edges) - 1));
            next_edge = viable_next_edges(next_edge_selection);

            % Check for scheduled stop
            if ismember(events(end,2), planned_stops(:,2)) && ~ismember(next_edge, planned_stops(:,2))
                if ismember(events(end,2), planned_stops(:,2)) && ismember(next_edge, planned_stops(:,2))
                    error("Planned stops must not be on adjacent edges.")
                end

                planned_stops(planned_stops(:,2) == events(end,2), :) = 0; % Remove the stop
                departure_timestep = min(edge_transition.timestep + params.dwell_timesteps, params.n_timesteps);
                [v_targets_working_set start_braking_timestep] = addBraking(params, global_speeds, v_targets_working_set, edge_transition, false, departure_timestep, 0);
            % Check for overspeed on entering new edge
            elseif abs(edge_transition.speed) > network.speed_limits(next_edge)
                [v_targets_working_set start_braking_timestep] = addBraking(params, global_speeds, v_targets_working_set, edge_transition, true, [], sign(edge_transition.speed) * network.speed_limits(next_edge));
            else 
                v_targets_modified = false;
            end
        end

        if v_targets_modified
            assert(size(v_targets, 2) == size(v_targets_working_set, 2));
            % Write new v targets for braking curve
            findVTargetIdxsOnEdge = @(x) find(x(:,1) >= start_braking_timestep & x(:,1) < edge_transition.timestep);

            v_targets(findVTargetIdxsOnEdge(v_targets), :) = [];
            v_targets = cat(1, v_targets, v_targets_working_set(findVTargetIdxsOnEdge(v_targets_working_set), :));

            assert(isequal(unique(v_targets(:,1), 'stable'), v_targets(:,1)));

            % Jump back to before earliest modified speed target point
            events = events(events(:, 1) < start_braking_timestep, :);
            if isempty(events)
                events(1,:) = [1, initial_position(1), initial_position(2), initial_position(3), initial_speed];
            end
            continue;
        end

        % If the transition needs no further modification finalize it

        % Write trajectory for this edge
        traj(events(end,1):edge_transition.timestep - 1, 1) = events(end, 2);
        traj(events(end,1):edge_transition.timestep - 1, 2) = edge_trajectory;
        traj(events(end,1):edge_transition.timestep - 1, 3) = events(end, 4);
        traj(events(end,1):edge_transition.timestep - 1, 4) = edge_speeds;

        % Write initial state for next edge
        edge_entrance_point = (network.edge_cols(next_edge) == edge_transition.traversed_node);
        new_event_entry_idx = size(events, 1) + 1;
        events(new_event_entry_idx, 1) = edge_transition.timestep;
        events(new_event_entry_idx, 2) = next_edge; % New Edge
        events(new_event_entry_idx, 3) = edge_entrance_point + (-1) * (edge_entrance_point*2 - 1) * abs(edge_transition.extra_movement / network.edge_values(next_edge)); % New Position on Edge
        % new_train_orientation      =         edge_entrance_point      XOR node_traversal_direction    ( XNOR is multiplication )
        events(new_event_entry_idx, 4) = (-1) * (edge_entrance_point * 2 - 1) * edge_transition.node_traversal_direction; % New Orientation on Edge
        events(new_event_entry_idx, 5) = edge_transition.speed;

        if abort > 1e2
            error("Trajectory construction failed.");
        end
        abort = abort + 1;
    end
end
