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

    events = assignEdgeTransitions(network, params, solution, initial_position, initial_speed);

    % TODO: write this function
end

function events = assignEdgeTransitions(network, params, solution, initial_position, initial_speed)
    %% This is the event-based part of the simulation
    % Also enforces dwell times for scheduled edges

    % events dimensions (edge_changes, 4)
    % events values (timestep, new_edge, position_on_edge, train_orientation_on_edge)

    acceleration = interpolateSolutionCurve(solution(1:params.n_acceleration_vars)*params.n_timesteps, solution(params.n_acceleration_vars+1:2*params.n_acceleration_vars), 1:params.n_timesteps);
    switch_directions = solution(2 * params.n_acceleration_vars + 1:2 * params.n_acceleration_vars + params.n_switch_vars);
    position = generatePositionCurve(params, acceleration, initial_speed);

    events(1, 1) = 1; 
    events(1, 2) = initial_position(1);
    events(1, 3) = initial_position(2);
    events(1, 4) = initial_position(3);

    pivot_timestep = 1;
    i_edge_change = 1;
    while pivot_timestep < params.n_timesteps
        % Determine position on edge
        current_edge_length = network.edge_values(events(i_edge_change, 2));
        remaining_backward_length = events(i_edge_change, 3) * current_edge_length;
        remaining_forward_length = current_edge_length - remaining_backward_length;
        base_position = position(pivot_timestep);

        % Find next edge exit 
        edge_trajectory = events(i_edge_change, 4) * (position(pivot_timestep:end) - base_position);
        next_pivot_timestep = (pivot_timestep - 1) + find((edge_trajectory > remaining_forward_length | edge_trajectory < -remaining_backward_length), 1);
        edge_exit_point = any((edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) > remaining_forward_length));
        if isempty(next_pivot_timestep)
            next_pivot_timestep = params.n_timesteps;
        end
        assert(next_pivot_timestep <= params.n_timesteps);

        % Leave current edge
        % node_traversal_direction = edge_exit_point XNOR old_train_orientation
        if edge_exit_point
            traversed_node = network.edge_cols(events(i_edge_change, 2));
            node_traversal_direction = events(i_edge_change, 4);
            extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) - remaining_forward_length;
        else
            traversed_node = network.edge_rows(events(i_edge_change, 2));
            node_traversal_direction = -events(i_edge_change, 4);
            extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) + remaining_backward_length;
        end

        viable_next_edges = network.adjacent_edge_list{traversed_node};
        viable_next_edges = viable_next_edges(viable_next_edges~=events(i_edge_change, 2));

        if isempty(viable_next_edges)
            % Stay stationary until trajectory comes back around
            % TODO: proper braking curve
            % Use the old base position and forw/backw lengths
            deadend_edge_trajectory = events(i_edge_change, 4) * (position(next_pivot_timestep:end) - base_position);
            if edge_exit_point
                deadend_next_pivot_timestep = (next_pivot_timestep - 1) + find(deadend_edge_trajectory < remaining_forward_length, 1);
            else
                deadend_next_pivot_timestep = (next_pivot_timestep - 1) + find(deadend_edge_trajectory > -remaining_backward_length, 1);
            end
            if isempty(deadend_next_pivot_timestep)
                deadend_next_pivot_timestep = params.n_timesteps;
            end
            assert(deadend_next_pivot_timestep <= params.n_timesteps);

            % Stay stationary at dead end of edge 
            events(i_edge_change + 1, 1) = deadend_next_pivot_timestep;
            events(i_edge_change + 1, 2) = events(i_edge_change, 2);
            events(i_edge_change + 1, 3) = double(edge_exit_point);
            events(i_edge_change + 1, 4) = events(i_edge_change, 4);

            pivot_timestep = deadend_next_pivot_timestep;
        else
            % Decide next edge
            next_edge_selection = 1 + round(switch_directions(i_edge_change) * (length(viable_next_edges) - 1));
            next_edge = viable_next_edges(next_edge_selection);

            % TODO: Check for a scheduled stop and adjust curve

            edge_entrance_point = (network.edge_cols(next_edge) == traversed_node);
            events(i_edge_change + 1, 1) = next_pivot_timestep;
            events(i_edge_change + 1, 2) = next_edge; % New Edge
            events(i_edge_change + 1, 3) = edge_entrance_point + (-1) * (edge_entrance_point*2 - 1) * abs(extra_movement / network.edge_values(next_edge)); % New Position on Edge
            % new_train_orientation      =         edge_entrance_point      XOR node_traversal_direction    ( XNOR is multiplication )
            events(i_edge_change + 1, 4) = (-1) * (edge_entrance_point*2 - 1) * node_traversal_direction; % New Orientation on Edge

            pivot_timestep = next_pivot_timestep;
        end

        i_edge_change = i_edge_change + 1;
    end
end

function traj = assignTrajectory()
    %% Combines position curve and edge transitions into a trajectory on graph

    % TODO: write this function

    % trajectory dimensions (3, timestep)
    % trajectory values (edge 0-n, position on edge 0-1, train orientation on edge -1, 1)

    traj(1, :) = initial_position(1);
    traj(2, :) = initial_position(2);
    traj(3, :) = initial_position(3);

    % Set trajectory on the current edge
    traj(1, pivot_timestep:next_pivot_timestep-1) = traj(1, pivot_timestep);
    traj(2, pivot_timestep:next_pivot_timestep-1) = traj(2, pivot_timestep) + edge_trajectory(1:next_pivot_timestep - (pivot_timestep-1) - 1) / current_edge_length;
    traj(3, pivot_timestep:next_pivot_timestep-1) = traj(3, pivot_timestep);

    traj(:, params.n_timesteps) = traj(:, params.n_timesteps-1);
end

function position = generatePositionCurve(params, acceleration, initial_speed)
    acceleration =  (acceleration * 2 - 1) * params.max_accel;
    speeds = initial_speed + cumsum(acceleration);
    acceleration(speeds>params.max_speed & acceleration>0) = 0;
    acceleration(speeds<-params.max_speed & acceleration<0) = 0;
    speeds = initial_speed + cumsum(acceleration);
    position = cumsum(speeds);
end

function y_new = interpolateSolutionCurve(x, y, x_new)
    %% Interpolate sparse curve representation to continuous one and normalize
    [~ , unique_idxs, ~] = unique(x);
    y_new = interp1(x(unique_idxs), y(unique_idxs), x_new, 'linear', 'extrap');
    y_new(y_new>1) = 1;
    y_new(y_new<0) = 0;
end