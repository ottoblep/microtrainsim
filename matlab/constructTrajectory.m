function [traj,arrival_events] = constructTrajectory(network, params, solution, initial_position, initial_speed, planned_stops)
    %% Constructs a single train trajectory on the graph

    % Acceleration vars are a sparse representation of the acceleration curves.
    % They are made up of support points (time-value pairs) that are then interpolated during trajectory construction.
    % This reduces the variable count and produces less erratic movement.
    % Switch directions are represent the next edge when hitting a node 
    % Planned stops are made automatically when the corresponding edge is passed and speed is slow enough

    % solution dimensions ( 2 * n_acceleration_vars + n_switch_vars )
    % solution values (accel_time 0-1, accel_value 0-1, switch_direction 0-1)
    % solution serialization (accel_timesteps, accel_values, sw_direction)

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

    acceleration = interpolateSolutionCurve(solution(1:params.n_acceleration_vars)*params.n_timesteps, solution(params.n_acceleration_vars+1:2*params.n_acceleration_vars), 1:params.n_timesteps);
    acceleration =  (acceleration * 2 - 1) * params.max_accel;
    switch_directions = solution(2 * params.n_acceleration_vars + 1:2 * params.n_acceleration_vars + params.n_switch_vars);
    speeds = initial_speed + cumsum(acceleration);
    speeds(speeds>params.max_speed) = params.max_speed;
    speeds(speeds<-params.max_speed) = -params.max_speed;
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
            [position, speeds, start_braking_timestep] = addStop(params, position, speeds, next_pivot_timestep, deadend_next_pivot_timestep, initial_speed, initial_position, 0);

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
            if ismember(next_edge, planned_stops(:,2))
                departure_timestep = min(next_pivot_timestep + params.dwell_timesteps, params.n_timesteps);
                [position, speeds, start_braking_timestep] = addStop(params, position, speeds, next_pivot_timestep, departure_timestep, initial_speed, initial_position, 1);
                planned_stops(planned_stops(:,2) == next_edge, :) = 0; % Remove planned stop
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

function [position, speeds, start_braking_timestep] = addStop(params, position, speeds, arrival_timestep, departure_time, initial_speed, initial_position, overshoot)
    %% Modifies position curve to stop around a certain position indicated by a timestep on the old curve

    approach_direction = sign(speeds(arrival_timestep));
    pre_stop_timesteps = 1:arrival_timestep;
    % Only consider time since last direction change
    first_approach_idx = find(sign(speeds(1:arrival_timestep)) ~= approach_direction, 1, 'last') + 1;
    if isempty(first_approach_idx)
        first_approach_idx = 1;
    end

    % Select timestep to start braking (will over or undershoot due to discretization)
    % Discrete Formula for distance covered under k braking steps
    % k*v(n) - a*k*(k+1)/2 
    % maximum (stop point) at (v-a/2)/a
    v = speeds(first_approach_idx:arrival_timestep);
    k = floor(v./params.max_accel);
    d_brake = k .* v - params.max_accel * k .* (k+1) * 0.5;
    d = abs(position(first_approach_idx:arrival_timestep) - position(arrival_timestep));
    start_braking_timestep = first_approach_idx + find(d < d_brake, 1, 'first') + (2 * overshoot - 1);
    if isempty(start_braking_timestep)
        start_braking_timestep = 1;
    end

    % Calculate braking curve
    if start_braking_timestep == 1
        n_full_braking_steps = floor(abs(initial_speed)/params.max_accel);
        a_last_step = abs(initial_speed) - n_full_braking_steps * params.max_accel;
    else
        n_full_braking_steps = floor(abs(speeds(start_braking_timestep-1))/params.max_accel);
        a_last_step = abs(speeds(start_braking_timestep-1)) - n_full_braking_steps * params.max_accel;
    end
    if all(n_full_braking_steps + 1 > params.n_timesteps - start_braking_timestep)
        warning("Not enough time to brake.");
        return;
    end
    end_braking_timestep = start_braking_timestep + n_full_braking_steps;

    % Adjust acceleration values
    acceleration = [speeds(1)-initial_speed diff(speeds)];
    acceleration(start_braking_timestep:end_braking_timestep - 1) = -approach_direction * params.max_accel;
    acceleration(end_braking_timestep) = -approach_direction * a_last_step;
    acceleration(end_braking_timestep+1:departure_time) = 0;

    % Recalculate positions
    if start_braking_timestep == 1
        speeds = initial_speed + cumsum(acceleration);
        position = cumsum(speeds);
    else
        speeds(start_braking_timestep:params.n_timesteps) = speeds(start_braking_timestep - 1) + cumsum(acceleration(start_braking_timestep:params.n_timesteps));
        position(start_braking_timestep:params.n_timesteps) = position(start_braking_timestep - 1) + cumsum(speeds(start_braking_timestep:params.n_timesteps));
    end

    assert((first_approach_idx <= start_braking_timestep) && (start_braking_timestep <= end_braking_timestep) && (end_braking_timestep <= departure_time));
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

function y_new = interpolateSolutionCurve(x, y, x_new)
    %% Interpolate sparse curve representation to continuous one and normalize
    [~ , unique_idxs, ~] = unique(x);
    y_new = interp1(x(unique_idxs), y(unique_idxs), x_new, 'linear', 'extrap');
    y_new(y_new>1) = 1;
    y_new(y_new<0) = 0;
end