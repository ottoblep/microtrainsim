function [traj, events, n_fullfilled_stops] = constructTrajectory(network, params, solution, initial_position, initial_speed, planned_stops)
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
    [v_target_timesteps , unique_idxs, ~] = unique(v_target_timesteps);
    v_target_values = solution(params.n_v_target_vars + 1:2 * params.n_v_target_vars);
    v_target_values = (2 * v_target_values(unique_idxs) - 1) * params.max_speed;
    v_targets = [v_target_timesteps v_target_values];

    % We move from edge to edge (events are train state at first timestep on new edge)
    % p(n) - p(n-1) = v(n-1)
    % v(n) - v(n-1) = a(n-1)

    events(1) = [1, initial_position(1), initial_position(2), initial_position(3), initial_speed];
    traj(1,:) = events(1, 2:5); 

    while True 
        if (events(end, 5) > network.speed_limits(traj(pivot_timestep, 1) || events(end, 5) > params.max_speed)
            error("Speed limit could not be satisfied.");
        end

        % Adjust copy of speed targets for current edge speed limit
        v_targets_working_set = min(v_targets(2), network.speed_limits(events(end, 2));

        % Find next edge exit
        [edge_transition edge_trajectory] = simulateEdge(network, params, events(end, :), v_targets_working_set);

        if isempty(edge_transition)
            % Finish simulation
            traj(events(end,1):params.n_timesteps, 1) = events(end, 2);
            traj(events(end,1):params.n_timesteps, 2) = edge_trajectory;
            traj(events(end,1):params.n_timesteps, 3) = events(end, 4);
            traj(events(end,1):params.n_timesteps, 4) = events(end, 5);
            break;
        end

        viable_next_edges = network.adjacent_edge_list{traversed_node};
        viable_next_edges = viable_next_edges(viable_next_edges~=initial_edge_state(1));

        revisit_events = true;

        % Check for dead end
        if isemtpy(viable_next_edges)
            % Stop train until speed target is in the other direction
            departure_timestep = edge_transition.timestep - 1 + find(sign(speeds_edge(edge_transition.timestep:params.n_timesteps)) ~= node_traversal_direction, 1, 'first');
            % TODO: adjust speed targets
        else
            % Decide next edge
            next_edge_selection = 1 + round(switch_directions(i_edge_change) * (length(viable_next_edges) - 1));
            next_edge = viable_next_edges(next_edge_selection);

            % Check for scheduled stop
            if ismember(initial_edge_state(1), planned_stops(:,2)) && ~ismember(next_edge, planned_stops(:,2))
                if ismember(initial_edge_state(1), planned_stops(:,2)) && ismember(next_edge, planned_stops(:,2))
                    error("Planned stops must not be on adjacent edges.")
                end

                planned_stops(planned_stops(:,2) == initial_edge_state(1), :) = 0; % Remove the stop
                departure_timestep = min(edge_transition.timestep + params.dwell_timesteps, params.n_timesteps);
                % TODO: adjust speed targets
            % Check for overspeed on entering new edge
            elseif edge_transition.speed > network.speed_limits(next_edge)
                % TODO: adjust speed targets
            else 
                revisit_events = false;
            end
        end

        if revisit_events
            % Jump back to before earliest modified speed target point
            events = events(events(:, 1) < start_braking_timestep, :);
            if isempty(sim_events)
                events(1) = [1, initial_position(1), initial_position(2), initial_position(3), initial_speed];
            end
            continue;
        end

        % Write trajectory for this edge
        traj(events(end,1):edge_transition.timestep - 1, 1) = events(end, 2);
        traj(events(end,1):edge_transition.timestep - 1, 2) = edge_trajectory;
        traj(events(end,1):edge_transition.timestep - 1, 3) = events(end, 4);
        traj(events(end,1):edge_transition.timestep - 1, 4) = events(end, 5);

        % Write adjusted v_targets for this edge
        % TODO

        % Write initial state for next edge
        edge_entrance_point = (network.edge_cols(next_edge) == edge_transition.traversed_node);
        events(end + 1, 1) = edge_transition.timestep;
        events(end + 1, 2) = next_edge; % New Edge
        events(end + 1, 3) = edge_entrance_point + (-1) * (edge_entrance_point*2 - 1) * abs(edge_transition.extra_movement / network.edge_values(next_edge)); % New Position on Edge
        % new_train_orientation      =         edge_entrance_point      XOR node_traversal_direction    ( XNOR is multiplication )
        events(end + 1, 4) = (-1) * (edge_entrance_point * 2 - 1) * edge_transition.node_traversal_direction; % New Orientation on Edge
        events(end + 1, 5) = edge_transition.speed;
    end
end
