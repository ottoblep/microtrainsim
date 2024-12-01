function [edge_transition edge_trajectory] = simulateEdge(network, params, initial_edge_state, v_targets) 
        %% Continuous simulation of movement on one edge
        % Does not check any constraints

        % initial_edge_state values (timestep [0,n], edge [1-n], position on edge [0,1], train orientation on edge {-1, 1}, speed [-Inf,Inf])

        % v_target dimensions (n, 2)
        % v_target values (timestep [0,n], speed [-Inf,Inf])

        % Data starts at the initial edge state
        v_targets_cont_edge = interp1(v_targets(:,1), v_targets(:,2), initial_edge_state(1):params.n_timesteps, 'previous', 'extrap');

        edge_speeds = zeros(1, params.n_timesteps - initial_edge_state(1) + 1);
        edge_trajectory = zeros(1, params.n_timesteps - initial_edge_state(1) + 1);
        edge_trajectory(1) = initial_edge_state(3);
        edge_length = network.edge_values(initial_edge_state(2));
        edge_length_divisor = 1 / edge_length;
;
        edge_speeds(1) = initial_edge_state(5);
        edge_transition = [];
        for i = 2:params.n_timesteps
            speed_disparity = v_targets_cont_edge(i-1) - edge_speeds(i-1);
            edge_speeds(i) = edge_speeds(i-1) + sign(speed_disparity) * min(abs(speed_disparity), params.max_accel);
            edge_trajectory(i) = edge_trajectory(i-1) + (initial_edge_state(4) * edge_speeds(i-1) * edge_length_divisor);

            if edge_trajectory(i) > 1
                edge_exit_point = 1;
                edge_transition.timestep = i;
                break;
            elseif edge_trajectory(i) < 0
                edge_exit_point = 0;
                edge_transition.timestep = i;
                break;
            end
        end

        % If no other transition exists the simulation is finished
        if isempty(edge_transition) return; end

        % Determine transition properties 
        % node_traversal_direction = edge_exit_point XNOR old_train_orientation
        if edge_exit_point % exit forwards
            edge_transition.traversed_node = network.edge_cols(initial_edge_state(2));
            edge_transition.node_traversal_direction = initial_edge_state(4);
            edge_transition.extra_movement = (edge_trajectory(edge_transition.timestep) - 1) * edge_length;
        else
            edge_transition.traversed_node = network.edge_rows(initial_edge_state(2));
            edge_transition.node_traversal_direction = -initial_edge_state(4);
            edge_transition.extra_movement = -edge_trajectory(edge_transition.timestep) * edge_length;
        end

        edge_transition.speed = edge_speeds(edge_transition.timestep);
        edge_speeds(edge_transition.timestep:end) = [];
        edge_trajectory(edge_transition.timestep:end) = [];
        % Return to global reference frame
        edge_transition.timestep = edge_transition.timestep + initial_edge_state(1) - 1;
end
