function [traj,arrival_events] = constructTrajectory(network, params, solution, initial_position, initial_speed, planned_stops)
    %% Constructs a single train trajectory on the graph

    % Target velocity vars are a sparse representation of the target velocity curves.
    % They are made up of support points (time-value pairs) that are then interpolated during trajectory construction.
    % This reduces the variable count and produces less erratic movement.
    % Switch directions are represent the next edge when hitting a node 
    % Planned stops are made automatically when the corresponding edge is passed and speed is slow enough

    % solution dimensions ( 2 * n_v_target_vars + n_switch_vars )
    % solution values (v_target_time 0-1, v_target_value 0-1, switch_direction 0-1)
    % solution serialization (v_target_timesteps, v_target_values, sw_direction)

    % initial position dimensions (3)
    % initial position values (edge 0-n, position on edge 0-1, train orientation on edge -1, 1)

    % initial speed dimensions (1)
    % initial speed values (speed: -max_speed-max_speed)

    % planned stops dimenstions (n, 3)
    % planned stops values (edge 0-n, arrival_time_fraction_of_timesteps 0-1)

    % sim_events dimensions (edge_changes, 4)
    % sim_events values (timestep, new_edge, position_on_edge, train_orientation_on_edge)

    % trajectory dimensions (3, timestep)
    % trajectory values (edge 0-n, position on edge 0-1, train orientation on edge -1, 1)

    [sim_events, position] = assignEdgeTransitions(network, params, solution, initial_position, initial_speed, planned_stops);

    traj = assignTrajectory(network, params, position, sim_events, initial_position);

    arrival_events = sim_events(:, 1:2);
end

function [sim_events, position] = assignEdgeTransitions(network, params, solution, initial_position, initial_speed, planned_stops)
    %% This is the event-based part of the simulation
    % Also enforces stopping for scheduled edges and dead ends

    switch_directions = solution(2 * params.n_v_target_vars + 1:2 * params.n_v_target_vars + params.n_switch_vars);
    speeds = constructMovement(params, solution, initial_speed);
    position = cumsum(speeds);

    sim_events(1,:) = [1 initial_position(1) initial_position(2) initial_position(3)];

    pivot_timestep = 1;
    i_edge_change = 1;
    while pivot_timestep < params.n_timesteps
        % Determine position on edge
        current_edge_length = network.edge_values(sim_events(i_edge_change, 2));
        remaining_backward_length = sim_events(i_edge_change, 3) * current_edge_length;
        remaining_forward_length = current_edge_length - remaining_backward_length;

        % Find next edge exit 
        edge_trajectory = sim_events(i_edge_change, 4) * (position(pivot_timestep:end) - position(pivot_timestep));
        next_pivot_timestep = (pivot_timestep - 1) + find((edge_trajectory > remaining_forward_length | edge_trajectory < -remaining_backward_length), 1);
        edge_exit_point = any((edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) > remaining_forward_length));
        if isempty(next_pivot_timestep)
            return;
        end

        % Leave current edge
        % node_traversal_direction = edge_exit_point XNOR old_train_orientation
        if edge_exit_point
            traversed_node = network.edge_cols(sim_events(i_edge_change, 2));
            node_traversal_direction = sim_events(i_edge_change, 4);
            extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) - remaining_forward_length;
        else
            traversed_node = network.edge_rows(sim_events(i_edge_change, 2));
            node_traversal_direction = -sim_events(i_edge_change, 4);
            extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) + remaining_backward_length;
        end

        viable_next_edges = network.adjacent_edge_list{traversed_node};
        viable_next_edges = viable_next_edges(viable_next_edges~=sim_events(i_edge_change, 2));

        if isempty(viable_next_edges)
            % Find time train would turn around 
            deadend_next_pivot_timestep = next_pivot_timestep - 1 + find(sign(speeds(next_pivot_timestep:params.n_timesteps)) ~= node_traversal_direction, 1, 'first');
            if isempty(deadend_next_pivot_timestep)
                deadend_next_pivot_timestep = params.n_timesteps;
            end

            % Modify curve so dead end is no longer hit
            assert(next_pivot_timestep < deadend_next_pivot_timestep);
            [solution, start_braking_timestep] = addStop(params, position, speeds, solution, next_pivot_timestep, deadend_next_pivot_timestep, false);
            speeds = constructMovement(params, solution, initial_speed);
            position = cumsum(speeds);

            % Revisit some past events that may have changed timing 
            sim_events = sim_events(sim_events(:, 1) < start_braking_timestep, :);
            if isempty(sim_events)
                sim_events(1,:) = [1 initial_position(1) initial_position(2) initial_position(3)];
                pivot_timestep = 1;
            else
                pivot_timestep = sim_events(end, 1);
            end
            i_edge_change = size(sim_events, 1);
        else
            % Decide next edge
            next_edge_selection = 1 + round(switch_directions(i_edge_change) * (length(viable_next_edges) - 1));
            next_edge = viable_next_edges(next_edge_selection);

            % Check for a scheduled stop that has not yet been visited
            if ismember(next_edge, planned_stops(:,1))
                planned_stops(planned_stops(:,2) == next_edge, :) = 0; % Remove planned stop
                departure_timestep = min(next_pivot_timestep + params.dwell_timesteps, params.n_timesteps);

                [solution, start_braking_timestep] = addStop(params, position, speeds, solution, next_pivot_timestep, departure_timestep, true);
                speeds = constructMovement(params, solution, initial_speed);
                position = cumsum(speeds);

                % Revisit some past events that may have changed timing 
                sim_events = sim_events(sim_events(:, 1) < start_braking_timestep, :);
                if isempty(sim_events)
                    sim_events(1,:) = [1 initial_position(1) initial_position(2) initial_position(3)];
                    pivot_timestep = 1;
                else
                    pivot_timestep = sim_events(end, 1);
                end
                i_edge_change = size(sim_events, 1);
            end

            edge_entrance_point = (network.edge_cols(next_edge) == traversed_node);
            sim_events(i_edge_change + 1, 1) = next_pivot_timestep;
            sim_events(i_edge_change + 1, 2) = next_edge; % New Edge
            sim_events(i_edge_change + 1, 3) = edge_entrance_point + (-1) * (edge_entrance_point*2 - 1) * abs(extra_movement / network.edge_values(next_edge)); % New Position on Edge
            % new_train_orientation      =         edge_entrance_point      XOR node_traversal_direction    ( XNOR is multiplication )
            sim_events(i_edge_change + 1, 4) = (-1) * (edge_entrance_point*2 - 1) * node_traversal_direction; % New Orientation on Edge

            pivot_timestep = next_pivot_timestep;
            i_edge_change = i_edge_change + 1;
        end
    end
end

function traj = assignTrajectory(network, params, position, sim_events, initial_position)
    %% Combines position curve and edge transitions into a continuous trajectory on graph

    traj(1, :) = initial_position(1);
    traj(2, :) = initial_position(2);
    traj(3, :) = initial_position(3);

    for i_transition = 1:size(sim_events, 1)
        if i_transition < size(sim_events, 1)
            edge_length = network.edge_values(sim_events(i_transition, 2));
            edge_trajectory = sim_events(i_transition, 4) * (position(sim_events(i_transition, 1):sim_events(i_transition + 1, 1) - 1) - position(sim_events(i_transition, 1)));
            traj(1,sim_events(i_transition, 1):sim_events(i_transition + 1, 1) - 1) = sim_events(i_transition, 2);
            traj(2,sim_events(i_transition, 1):sim_events(i_transition + 1, 1) - 1) = sim_events(i_transition, 3) + edge_trajectory / edge_length;
            traj(3,sim_events(i_transition, 1):sim_events(i_transition + 1, 1) - 1) = sim_events(i_transition, 4);
        else
            edge_length = network.edge_values(sim_events(i_transition, 2));
            edge_trajectory = sim_events(i_transition, 4) * (position(sim_events(i_transition, 1):params.n_timesteps) - position(sim_events(i_transition, 1)));
            traj(1,sim_events(i_transition, 1):params.n_timesteps) = sim_events(i_transition, 2);
            traj(2,sim_events(i_transition, 1):params.n_timesteps) = sim_events(i_transition, 3) + edge_trajectory / edge_length;
            traj(3,sim_events(i_transition, 1):params.n_timesteps) = sim_events(i_transition, 4);
        end
    end
end

function [solution, start_braking_timestep] = addStop(params, position, speeds, solution, arrival_timestep, departure_time, overshoot)
    %% Modifies solution to stop around a certain position defined by a timestep on the old position curve

    approach_direction = sign(speeds(arrival_timestep));
    % Only consider time since last direction change
    first_approach_idx = find(sign(speeds(1:arrival_timestep-1)) ~= approach_direction, 1, 'last') + 1;
    if isempty(first_approach_idx)
        first_approach_idx = 1;
    end

    % Select timestep to start braking (will over or undershoot due to discretization)
    % Discrete Formula for distance covered under n full braking steps and one remainder braking step
    possible_braking_timesteps = first_approach_idx:arrival_timestep;
    distance_from_stop = abs(position(possible_braking_timesteps) - position(arrival_timestep));
    n = floor(abs(speeds(possible_braking_timesteps)) / params.max_accel);
    distance_covered_under_braking = n .* (abs(speeds(possible_braking_timesteps)) - 0.5 * (n+1) .* params.max_accel);
    position_error = distance_covered_under_braking - distance_from_stop;
    
    if overshoot
        subset_braking_idxs = find(position_error > 0);
    else
        subset_braking_idxs = find(position_error < 0);
    end
    if isempty(subset_braking_idxs)
        warning("Failed to overshoot on braking.");
        subset_braking_idxs = 1:numel(position_error);
    end

    [~, best_subset_braking_idx] = min(abs(position_error(subset_braking_idxs)));
    start_braking_timestep = first_approach_idx +  subset_braking_idxs(best_subset_braking_idx);

    if isempty(start_braking_timestep)
        start_braking_timestep = first_approach_idx;
    end
    
    % Adjust target velocity points
    v_target_timesteps = solution(1:params.n_v_target_vars);
    v_target_values = solution(params.n_v_target_vars + 1:2 * params.n_v_target_vars);
    [v_target_timesteps, v_target_sorted_idxs] = sort(v_target_timesteps);
    v_target_values = v_target_values(v_target_sorted_idxs);

    % Stay stationary 
    idxs_v_targets_during_stop = (v_target_timesteps >= start_braking_timestep) & (v_target_timesteps <= departure_time);
    v_target_values(idxs_v_targets_during_stop) = 0.5;
    % Start braking instantly
    idx_v_target_to_shift = find(v_target_timesteps > start_braking_timestep, 1, 'first');
    v_target_timesteps(idx_v_target_to_shift) = start_braking_timestep;
    v_target_values(idx_v_target_to_shift) = 0.5;

    solution(1:params.n_v_target_vars) = v_target_timesteps;
    solution(params.n_v_target_vars + 1:2 * params.n_v_target_vars) = v_target_values;
end

function speeds = constructMovement(params, solution, initial_speed)
    %% Constructs physically possible speed and position curves from target speeds points
    v_target_timesteps = solution(1:params.n_v_target_vars) * params.n_timesteps;
    [~, first_v_target_timestep_idx] = min(v_target_timesteps);
    v_target_timesteps(first_v_target_timestep_idx) = 1;
    [v_target_timesteps , unique_idxs, ~] = unique(v_target_timesteps);

    v_target_values = solution(params.n_v_target_vars + 1:2 * params.n_v_target_vars);
    v_target_values = (2 * v_target_values(unique_idxs) - 1) * params.max_speed;
    
    v_target = interp1(v_target_timesteps, v_target_values, 1:params.n_timesteps, 'previous');

    speeds = zeros(1, params.n_timesteps);
    speeds(1) = initial_speed;
    for i = 2:params.n_timesteps
        diff = v_target(i) - speeds(i-1);
        speeds(i) = speeds(i-1) + sign(diff) * min(abs(diff), params.max_accel);
    end

    % clf; hold on;
    % scatter(v_target_timesteps, v_target_values);
    % plot(v_target);
    % plot(speeds);
end