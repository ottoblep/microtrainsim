function [traj, events] = constructTrajectory(network, params, solution, initial_position, initial_speed)
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
    % trajectory dimensions (3, timestep)
    % trajectory values (edge 0-n, position on edge 0-1, train orientation on edge -1, 1)
    % events dimensions (2, n)
    % events values (edge, edge_entrance_timestep)

    acceleration = interpolateSolutionCurve(solution(1:n_acceleration_vars)*params.n_timesteps, solution(n_acceleration_vars+1:2*n_acceleration_vars), 1:params.n_timesteps);
    speeds = initial_speed + cumtrapz(((acceleration * 2) - 1) * max_accel);
    speeds(speeds>max_speed) = max_speed;
    speeds(speeds<-max_speed) = -max_speed;
    position = cumtrapz(speeds);

    switch_directions = solution(2 * n_acceleration_vars + 1:2 * n_acceleration_vars + n_switch_vars);

    events = [];
    traj(1, :) = initial_position(1);
    traj(2, :) = initial_position(2);
    traj(3, :) = initial_position(3);

    % Record possible initial placement on node
    if initial_position(2) == 0
        events(:, end+1) = [network.edge_rows(initial_position(1)), 1];
    elseif initial_position(2) == 1
        events(:, end+1) = [network.edge_cols(initial_position(1)), 1];
    end

    pivot_timestep = 1;
    i_edge_change = 1;
    while pivot_timestep < params.n_timesteps
        current_edge_length = network.edge_values(traj(1, pivot_timestep));
        remaining_backward_length = traj(2, pivot_timestep) * current_edge_length;
        remaining_forward_length = current_edge_length - remaining_backward_length;
        base_position = position(pivot_timestep);

        % Find next edge change
        edge_trajectory = traj(3, pivot_timestep) * (position(pivot_timestep:end) - base_position);
        next_pivot_timestep = (pivot_timestep - 1) + find((edge_trajectory > remaining_forward_length | edge_trajectory < -remaining_backward_length), 1);
        edge_exit_point = any((edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) > remaining_forward_length));
        if isempty(next_pivot_timestep)
            next_pivot_timestep = params.n_timesteps;
        end
        assert(next_pivot_timestep <= params.n_timesteps);

        % Set trajectory on the current edge
        traj(1, pivot_timestep:next_pivot_timestep-1) = traj(1, pivot_timestep);
        traj(2, pivot_timestep:next_pivot_timestep-1) = traj(2, pivot_timestep) + edge_trajectory(1:next_pivot_timestep - (pivot_timestep-1) - 1) / current_edge_length;
        traj(3, pivot_timestep:next_pivot_timestep-1) = traj(3, pivot_timestep);

        % Leave current edge
        % node_traversal_direction = edge_exit_point XNOR old_train_orientation
        if edge_exit_point
            traversed_node = network.edge_cols(traj(1, pivot_timestep));
            node_traversal_direction = traj(3, pivot_timestep);
            extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) - remaining_forward_length;
        else
            traversed_node = network.edge_rows(traj(1, pivot_timestep));
            node_traversal_direction = -traj(3, pivot_timestep);
            extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) + remaining_backward_length;
        end

        events(:, end+1) = [traversed_node, next_pivot_timestep];

        % Decide on next edge
        viable_next_edges = network.adjacent_edge_list{traversed_node};
        viable_next_edges = viable_next_edges(viable_next_edges~=traj(1, pivot_timestep));

        if isempty(viable_next_edges)
            % Break and stay stationary until trajectory comes back

        else
            next_edge_selection = 1 + round(switch_directions(i_edge_change) * (length(viable_next_edges) - 1));
            i_edge_change = i_edge_change + 1;
            next_edge = viable_next_edges(next_edge_selection);

            % Check for a scheduled stop

            edge_entrance_point = (network.edge_cols(next_edge) == traversed_node);
            % New Edge
            traj(1, next_pivot_timestep) = next_edge;
            % New Position on Edge
            traj(2, next_pivot_timestep) = edge_entrance_point + (-1) * (edge_entrance_point*2 - 1) * abs(extra_movement / network.edge_values(next_edge));
            % New Orientation on Edge
            % new_train_orientation      =         edge_entrance_point      XOR node_traversal_direction    ( XNOR is multiplication )
            traj(3, next_pivot_timestep) = (-1) * (edge_entrance_point*2 - 1) * node_traversal_direction;

            pivot_timestep = next_pivot_timestep;
        end
    end

    traj(:, params.n_timesteps) = traj(:, params.n_timesteps-1);
end

function y_new = interpolateSolutionCurve(x, y, x_new)
    %% Interpolate sparse curve representation to continuous one and normalize
    [~ , unique_idxs, ~] = unique(x);
    y_new = interp1(x(unique_idxs), y(unique_idxs), x_new, 'linear', 'extrap');
    y_new(y_new>1) = 1;
    y_new(y_new<0) = 0;
end