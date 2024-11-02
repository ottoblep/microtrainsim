function distance = trainDistance(network, traj_set, i_train, j_train, timestep)
    %% Calculates distance of two trains given a set of trajectories
    % trajectory_set dimensions (n_trains, 3, timestep)
    % trajectory_set values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)

    edge_i = int32(traj_set(i_train, 1, timestep));
    edge_j = int32(traj_set(j_train, 1, timestep));

    i_edge_length = network.edge_values(edge_i);
    j_edge_length = network.edge_values(edge_j);

    if edge_i == edge_j
        distance = abs(traj_set(i_train, 2, timestep) - traj_set(j_train, 2, timestep)) * i_edge_length;
    else
        % Find shortest path with precomputed distance matrix
        i_node_backward = network.edge_rows(edge_i);
        i_node_forward = network.edge_cols(edge_i);
        j_node_backward = network.edge_rows(edge_j);
        j_node_forward = network.edge_cols(edge_j);

        i_remaining_backward_length = traj_set(i_train, 2, timestep) * i_edge_length;
        i_remaining_forward_length = i_edge_length - i_remaining_backward_length;
        j_remaining_backward_length = traj_set(j_train, 2, timestep) * j_edge_length;
        j_remaining_forward_length = j_edge_length - j_remaining_backward_length;

        dist1 = i_remaining_backward_length + network.all_shortest_paths(i_node_backward,j_node_backward) + j_remaining_backward_length;
        dist2 = i_remaining_forward_length + network.all_shortest_paths(i_node_forward,j_node_backward) + j_remaining_backward_length;
        dist3 = i_remaining_backward_length + network.all_shortest_paths(i_node_backward,j_node_forward) + j_remaining_forward_length;
        dist4 = i_remaining_forward_length + network.all_shortest_paths(i_node_forward,j_node_forward) + j_remaining_forward_length;

        distance = min([dist1, dist2, dist3, dist4]);
    end
end