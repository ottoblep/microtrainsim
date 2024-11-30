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

    % We move from edge-change to edge-change (identified by first timestep on new edge)
    % p(n) - p(n-1) = v(n-1)

    events(1) = [1, initial_position(1), initial_position(2), initial_position(3), initial_speed];
    traj(1,:) = events(1, 2:5); 
    while True 
        if (events(end, 5) > network.speed_limits(traj(pivot_timestep, 1) || events(end, 5) > params.max_speed)
            error("Speed limit could not be satisfied.");
        end

        % Adjust copy of speed targets for current speed limit
        v_targets_working_set = v_targets();

        % Find next edge exit
        edge_transition = identifyNextEdgeExit(network, current_edge_state, v_targets, start_timestep, max_accel);
        if isempty(edge_transition)
            % TODO
        end

        viable_next_edges = network.adjacent_edge_list{traversed_node};
        viable_next_edges = viable_next_edges(viable_next_edges~=current_edge_state(1));

        revisit_events = true;

        % Check for dead end
        if isemtpy(viable_next_edges)
            % Stop train until speed target is in the other direction
            departure_timestep = next_pivot_timestep - 1 + find(sign(speeds(next_pivot_timestep:params.n_timesteps)) ~= node_traversal_direction, 1, 'first');
            % TODO: adjust speed targets
        else
            % Decide next edge
            next_edge_selection = 1 + round(switch_directions(i_edge_change) * (length(viable_next_edges) - 1));
            next_edge = viable_next_edges(next_edge_selection);

            % Check for scheduled stop
            if ismember(current_edge_state(1), planned_stops(:,2)) && ~ismember(next_edge, planned_stops(:,2))
                if ismember(current_edge_state(1), planned_stops(:,2)) && ismember(next_edge, planned_stops(:,2))
                    error("Planned stops must not be on adjacent edges.")
                end

                planned_stops(planned_stops(:,2) == current_edge_state(1), :) = 0; % Remove the stop
                departure_timestep = min(next_pivot_timestep + params.dwell_timesteps, params.n_timesteps);
                % TODO: adjust speed targets
            % Check for overspeed on entering new edge
            else if edge_transition(5) > network.speed_limits(next_edge)
                % TODO: adjust speed targets
            else 
                revisit_events = false;
            end
        end

        if revisit_events
            % Jump back to before earliest modified speed target point
            % TODO
        end

        % Write new event and trajectory
        % TODO
    end
end

function edge_transition = identifyNextEdgeExit(network, current_edge_state, v_targets, start_timestep, max_accel) 
        % Measure current edge
        current_edge_length = network.edge_values(current_edge_state(1));
        remaining_backward_length = current_edge_state(2) * current_edge_length;
        remaining_forward_length = current_edge_length - remaining_backward_length;

        % Construct preliminary trajectory on edge
        speeds = constructMovement(v_targets, start_timestep, params.n_timesteps, current_edge_state(4), params.max_accel);
        edge_trajectory = current_edge_state(3) * cumsum(speeds);

        % Find next edge exit 
        new_edge_timestep = find((edge_trajectory > remaining_forward_length | edge_trajectory < -remaining_backward_length), 1);
        edge_exit_point = any((edge_trajectory(new_edge_timestep - (pivot_timestep - 1)) > remaining_forward_length));

        % If no other transition exists the simulation is finished
        if isempty(new_edge_timestep)
            edge_transition = [];
            return;
        end

        % Determine transition properties 
        % node_traversal_direction = edge_exit_point XNOR old_train_orientation
        if edge_exit_point % exit forwards
            traversed_node = network.edge_cols(current_edge_state(1));
            node_traversal_direction = current_edge_state(3);
            extra_movement = edge_trajectory(new_edge_timestep) - remaining_forward_length;
        else
            traversed_node = network.edge_rows(current_edge_state(1));
            node_traversal_direction = -current_edge_state(3);
            extra_movement = edge_trajectory(new_edge_timestep) + remaining_backward_length;
        end

        edge_transition = [new_edge_timestep, traversed_node, node_traversal_direction, extra_movement, speeds(new_edge_timestep)];
end

function speeds = constructMovement(v_targets, start_timestep, end_timestep, initial_speed, max_accel)
    %% Constructs physically possible speed curve from target speeds points

    v_target = interp1(v_target_timesteps, v_target_values, start_timestep:end_timestep, 'previous', 'extrap');
    v_target(1) = v_target(1) - (start_timestep - 1);

    speeds = zeros(1, end_timestep - start_timestep + 1);
    speeds(1) = inital_speed;
    for i = 2:(end_timestep - start_timestep + 1)
        disparity = v_target(i) - speeds(i-1);
        speeds(i) = speeds(i-1) + sign(disparity) * min(abs(disparity), max_accel);
    end

    % clf; hold on;
    % scatter(v_target_timesteps, v_target_values);
    % plot(v_target);
    % plot(speeds);
end


% --------------------------------------------------------- OLD CODE

        if revisit_events
            % Reevaluate events from where position curve was modified
            sim_events = sim_events(sim_events(:, 1) < start_braking_timestep, :);
            if isempty(sim_events)
                sim_events(1,:) = [1 initial_position(1) initial_position(2) initial_position(3)];
                pivot_timestep = 1;
            else
                pivot_timestep = sim_events(end, 1);
            end
            i_edge_change = size(sim_events, 1);
            continue;
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

function [position, speeds, start_braking_timestep] = addStop(params, position, speeds, solution, arrival_timestep, departure_time, initial_speed)
    %% Modifies position curve to stop around a certain position defined by a timestep on the old position curve

    start_braking_timestep = findBrakingTimestep(position, speeds, arrival_timestep, params.max_accel);
    [v_target_timesteps, v_target_values] = extractSpeedTargetPoints(params, solution);

    % Consider only points after start of braking
    v_target_relevant_idxs = find(v_target_timesteps >= start_braking_timestep);
    v_target_values = v_target_values(v_target_relevant_idxs);
    v_target_timesteps = v_target_timesteps(v_target_relevant_idxs);

    % Stay stationary during dwell time
    idxs_v_targets_during_stop = (v_target_timesteps >= start_braking_timestep) & (v_target_timesteps <= departure_time);
    v_target_values(idxs_v_targets_during_stop) = 0;

    % Place speed target point of 0 at the start of braking 
    v_target_values(v_target_timesteps == start_braking_timestep) = [];
    v_target_timesteps(v_target_timesteps == start_braking_timestep) = [];
    v_target_timesteps(end+1) = start_braking_timestep;
    v_target_values(end+1) = 0;

    % Recalculate position curve
    if numel(v_target_timesteps) == 1
        v_target = ones(1, params.n_timesteps - start_braking_timestep + 1) * v_target_values(1);
    else
        v_target = interp1(v_target_timesteps, v_target_values, start_braking_timestep:params.n_timesteps, 'previous', 'extrap');
    end
    
    % clf; hold on;
    % plot(position,'DisplayName',"position old");
    % plot(speeds,'DisplayName',"speeds old");

    for i = start_braking_timestep:params.n_timesteps
        if i == 1
            disparity = v_target(i - start_braking_timestep + 1) - initial_speed;
            speeds(i) = initial_speed + sign(disparity) * min(abs(disparity), params.max_accel);
        else
            disparity = v_target(i - start_braking_timestep + 1) - speeds(i-1);
            speeds(i) = speeds(i-1) + sign(disparity) * min(abs(disparity), params.max_accel);
        end
    end

    if start_braking_timestep == 1
        position(start_braking_timestep:params.n_timesteps) = cumsum(speeds(start_braking_timestep:params.n_timesteps));
    else
        position(start_braking_timestep:params.n_timesteps) = position(start_braking_timestep - 1) + cumsum(speeds(start_braking_timestep:params.n_timesteps));
    end

    % plot(position,'DisplayName', "position new");
    % plot(speeds,'DisplayName', "speeds new");
    % scatter(start_braking_timestep, 400,'DisplayName', "start braking timestep");
    % scatter(v_target_timesteps, v_target_values, 'DisplayName',"new speed targets");
    % legend();
end

function start_braking_timestep = findBrakingTimestep(position, speeds, arrival_timestep, max_accel)
    %% Calculate start of braking in order to stop at a position defined by a timestep on the position curve

    approach_direction = sign(speeds(arrival_timestep));
    % Only consider time since last direction change
    first_approach_idx = find(sign(speeds(1:arrival_timestep-1)) == -approach_direction, 1, 'last') + 1;
    if isempty(first_approach_idx)
        first_approach_idx = 1;
    end

    % Select timestep to start braking (will purposefully undershoot due to discretization)
    % Discrete Formula for distance covered under n full braking steps and one remainder braking step
    possible_braking_timesteps = first_approach_idx:arrival_timestep;
    distance_from_stop = abs(position(possible_braking_timesteps) - position(arrival_timestep));
    n = floor(abs(speeds(possible_braking_timesteps)) / max_accel);
    distance_covered_under_braking = n .* (abs(speeds(possible_braking_timesteps)) - 0.5 * (n+1) .* max_accel);
    position_error = distance_covered_under_braking - distance_from_stop;
    
    subset_braking_idxs = find(position_error < 0);

    if isempty(subset_braking_idxs)
        %warning("Failed to undershoot on braking.");
        subset_braking_idxs = 1:numel(position_error);
    end

    [~, best_subset_braking_idx] = min(abs(position_error(subset_braking_idxs)));
    start_braking_timestep = first_approach_idx - 1 + subset_braking_idxs(best_subset_braking_idx);

    if isempty(start_braking_timestep)
        start_braking_timestep = first_approach_idx;
    end
end
