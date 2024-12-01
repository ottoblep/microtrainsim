function [network, params] = generateEnvironment(network_template)
    %% Network
    network.adjacency_matrix = readmatrix(strcat("./network_templates/", network_template, ".csv"));
    [network.edge_rows, network.edge_cols, network.edge_values] = find(network.adjacency_matrix);
    network.adjacent_edge_list = {};
    for node = 1:size(network.adjacency_matrix,1)
        network.adjacent_edge_list{node} = find((network.edge_rows == node) | (network.edge_cols == node));
    end
    network.speed_limits = readmatrix(strcat("./network_templates/", network_template, "_speed_limits.csv"));
    % All shortest path pairs
    tmp_adj = network.adjacency_matrix;
    tmp_adj(tmp_adj==0) = Inf;
    network.all_shortest_paths = distances(graph(tmp_adj,'upper'), 'Method', 'positive');
        

    %% Simulation Parameters
    params.n_timesteps = 200; % 10s timesteps
    params.n_v_target_vars = round(params.n_timesteps / 10); % Support points for target velocity
    params.min_separation = 100; % m
    params.max_speed = 83; % m/10s = 20km/h
    params.max_accel = 46.27; % m/(10s)² = 0-100kmh in 1m
    params.max_changeover_time = 1440; % 4hrs
    params.train_capacity = 400; % 400 Passengers
    params.dwell_timesteps = 12; % 2 minutes
    params.n_switch_vars = params.max_speed * params.n_timesteps / min(tmp_adj, [], 'all'); % Bounded by max possible edge changes

    %% Train Parameters
    params.initial_positions = readmatrix(strcat("./network_templates/", network_template, "_initial_positions.csv"));
    params.n_trains = size(params.initial_positions, 1);
    params.initial_speeds = zeros(params.n_trains, 1);

    params.destinations = readmatrix(strcat("./network_templates/", network_template, "_destinations.csv"));
    params.planned_stops = readmatrix(strcat("./network_templates/", network_template, "_planned_stops.csv"));
end

