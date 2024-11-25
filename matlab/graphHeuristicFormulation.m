function greedyHeuristicFormulation()
    [network, params] = generateEnvironment("crossover");
    assert(all(network.edge_values > params.max_speed)); % Requirement for trajectory construction

    [solution, traj_set] = greedySearch(network, params, 15, 0.1);
    %[solution, traj_set] = geneticGlobalSearch(network, params);
    %[solution, falsetraj_set] = particleSwarmSearch(network, params);

    %[solution, traj_set] = repairHeuristic(network, params, solution, traj_set);
    
    [solution, traj_set] = refineSolution(network, params, solution, traj_set);

    final_collision_score = collisionPenalties(network, traj_set, params.min_separation, params.max_speed)
    csvwrite("network.csv", network.adjacency_matrix);
    csvwrite("trajectories_edges.csv", squeeze(traj_set(:,1,:)));
    csvwrite("trajectories_positions.csv", squeeze(traj_set(:,2,:)));
end

function [network, params] = generateEnvironment(network_template)
    %% Network
    network.adjacency_matrix = readmatrix(strcat("../network_templates/", network_template, ".csv"));
    network.station_edges = readmatrix(strcat("../network_templates/", network_template, "_stations.csv"));
    [network.edge_rows, network.edge_cols, network.edge_values] = find(network.adjacency_matrix);
    network.adjacent_edge_list = {};
    for node = 1:size(network.adjacency_matrix,1)
        network.adjacent_edge_list{node} = find((network.edge_rows == node) | (network.edge_cols == node));
    end
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
    params.initial_positions = readmatrix(strcat("../network_templates/", network_template, "_initial_positions.csv"));
    params.n_trains = size(params.initial_positions, 1);
    params.initial_speeds = zeros(params.n_trains, 1);

    params.destinations = readmatrix(strcat("../network_templates/", network_template, "_destinations.csv"));
    params.planned_stops = readmatrix(strcat("../network_templates/", network_template, "_planned_stops.csv"));
end

%% Solution Construction

function [traj_set, event_set, n_fullfilled_stops] = constructTrajectorySet(network, params, solution)
    %% Constructs a set of train trajectories on the graph

    % solution dimensions (n_trains, 2 * n_v_target_vars + n_switch_vars)

    % initial position dimensions (n_trains, 3)
    % initial speed dimensions (n_trains)

    % planned stops dimenstions (n, 3)
    % planned stops values (train, edge, arrival_time_fraction_of_timesteps)

    % traj_set dimensions (n_trains, 3, timestep)

    % event_set dimensions (n, 3))
    % event_set values (train, timestep, presence_at_node)

    n_trains = size(solution,1);
    event_set = [];
    traj_set = zeros(n_trains, 3, params.n_timesteps);
    n_fullfilled_stops = 0;
    for i_train = 1:n_trains
        [traj_set(i_train, :, :), new_events, n_fullfilled_stops_train] = constructTrajectory(network, params, solution(i_train,:), params.initial_positions(i_train, :), params.initial_speeds(i_train), params.planned_stops(params.planned_stops(:,1)==i_train, 2:3));
        new_events(:, 1) = i_train;
        event_set = cat(1, event_set, new_events);
        n_fullfilled_stops = n_fullfilled_stops + n_fullfilled_stops_train;
    end
end

%% Helper Functions

function [edge_idx edge_pos] = nodeToEdgePos(network, node_idx)
    out_edges = find(network.edge_rows == node_idx);
    if ~isempty(out_edges)
        edge_idx = out_edges(1);
        edge_pos = 0;
    else
        in_edges = find(network.edge_cols == node_idx);
        edge_idx = in_edges(1);
        edge_pos = 1;
    end
end

function adj = randomPlanarGraph(n)
    % adj = spalloc(n, n, 6*n*n); % for finite planar graphs the average degree is strictly less than 6
    del = delaunay(randn(n, 2));
    for i_triangle = 1:size(del, 1)
        adj(del(i_triangle, 1),del(i_triangle, 2)) = 1;
        adj(del(i_triangle, 2),del(i_triangle, 3)) = 1;
        adj(del(i_triangle, 3),del(i_triangle, 1)) = 1;
    end
end

%% Objective Evalutation

function collision_score = collisionPenalties(network, traj_set, min_separation, max_speed)
    %% Evaluates a set of train trajectories

    % trajectory_set dimensions (n_trains, 3, timestep)
    % trajectory_set values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)

    collision_score = 0;

    % Separation Penalties
    % For each train pair update the minimum time to collision then skip that time and check again
    n_trains = size(traj_set, 1);
    n_timesteps = size(traj_set, 3);
    
    for i_train = 1:n_trains
        for j_train = i_train+1:n_trains
            timestep = 1;
            while timestep < n_timesteps
                distance = trainDistance(network, traj_set, i_train, j_train, timestep);

                if distance > min_separation
                    guaranteed_safe_time = int32(floor((distance - min_separation) / (2 * max_speed))) + 1;
                else
                    % Exponential penalty for closeness beyond the minimum separation
                    collision_score = collision_score - min(1e8, (min_separation/distance - 1));
                    guaranteed_safe_time = 1;
                end

                timestep = timestep + guaranteed_safe_time;
            end
        end
    end
end

function destination_score = destinationPenalties(network, traj_set, destinations)
    %% Penalizes train's distance from destination and speed>0 at the end of the trajectory

    % destinations dimensions (1, n_nodes)
    % destinations values (destination_node_idx)

    destination_score = 0;

    for i_train = 1:size(traj_set, 1)
        edge_i = int32(traj_set(i_train, 1, size(traj_set, 3)));
        i_edge_length = network.edge_values(edge_i);
        % Find shortest path with precomputed distance matrix
        i_node_backward = network.edge_rows(edge_i);
        i_node_forward = network.edge_cols(edge_i);
        i_remaining_backward_length = traj_set(i_train, 2, size(traj_set, 3)) * i_edge_length;
        i_remaining_forward_length = i_edge_length - i_remaining_backward_length;
        dist1 = i_remaining_backward_length + network.all_shortest_paths(i_node_backward, destinations(i_train));
        dist2 = i_remaining_forward_length + network.all_shortest_paths(i_node_forward, destinations(i_train));
        distance = min([dist1, dist2]);
        destination_score = destination_score - distance / 10;
    end
end

function [demand_score, transfer_graph_digraph, edge_flows] = demandSatisfaction(network, event_set, demand_matrix, max_changeover_time, train_capacity)
    %% Evaluates satisfied node-to-node demand for a given set of arrival/departure timesteps 

    % event_set dimenstions (3, n))
    % event_set values (presence_at_node, timestep, train)

    assert(all(size(demand_matrix) == size(network.adjacency_matrix)));

    if isempty(event_set)
        demand_score = 0;
        edge_flows = zeros(numel(find(network.adjacency_matrix)), 1);
        transfer_graph_digraph = digraph([]);
        return;
    end

    transfer_graph = constructTransferGraph(network, event_set, max_changeover_time, train_capacity);
    transfer_graph_digraph = digraph(transfer_graph);
    
    [flow_value, edge_flows] = maxMulticommodityFlowApprox(transfer_graph, transfer_graph_digraph, size(demand_matrix,1), demand_matrix, 0.2);    
    %[flow_value, edge_flows] = maxMulticommodityFlowLP(transfer_graph,
    %transfer_graph_digraph, size(demand_matrix,1), demand_matrix);

    demand_score = flow_value / sum(demand_matrix,'all');
end

function combined_score = combinedObjective(network, params, solution)
    [traj_set, event_set, n_fullfilled_stops] = constructTrajectorySet(network, params, solution);
    combined_score = collisionPenalties(network, traj_set, params.min_separation, params.max_speed);
    combined_score = combined_score + destinationPenalties(network, traj_set, params.destinations);
    combined_score = combined_score + n_fullfilled_stops * 10;
end

%% Solution Generation 

function solution = randomSolution(params)
    solution = rand(params.n_trains, 2 * params.n_v_target_vars + params.n_switch_vars);
end

function [solution, traj_set, round_best_score] = greedySolution(network, params, per_train_stall_time)
    %% Place each train greedily considering whole objective
    traj_set = [];
    event_set = [];
    solution = [];
    n_fullfilled_stops = 0;

    for i_train = 1:params.n_trains
        stall_timer = tic;
        round_best_traj_set = [];
        round_best_event_set = [];
        round_best_train_solution = [];
        round_best_score = -Inf;

        while toc(stall_timer) < per_train_stall_time
            new_train_solution = rand(1, 2 * params.n_v_target_vars + params.n_switch_vars);
            
            % Generate new trajectory
            [new_traj, new_events, n_fullfilled_stops_train] = constructTrajectory(network, params, new_train_solution, params.initial_positions(i_train, :), params.initial_speeds(i_train), params.planned_stops(params.planned_stops(:,1)==i_train, 2:3));

            % Add new trajectory to old set
            new_traj_set = cat(1, traj_set, reshape(new_traj, 1, size(new_traj,1), size(new_traj,2)));
            new_events(:, 1) = i_train;
            new_event_set = cat(1, event_set, new_events);

            % Test new trajectory set
            score = collisionPenalties(network, new_traj_set, params.min_separation, params.max_speed);
            score = score + destinationPenalties(network, new_traj_set, params.destinations);
            score = score + n_fullfilled_stops + n_fullfilled_stops_train * 10;

            if score > round_best_score
                round_best_score = score;
                stall_timer = tic;
                round_best_train_solution = new_train_solution;
                round_best_traj_set = new_traj_set;
                round_best_event_set = new_event_set;
            end
        end

        n_fullfilled_stops = n_fullfilled_stops + n_fullfilled_stops_train;
        traj_set = round_best_traj_set;
        event_set = round_best_event_set;
        solution = cat(1, solution, round_best_train_solution);
    end
end

%% Search Methods

function [solution, traj_set] = greedySearch(network, params, overall_stall_time, solution_stall_time)
    %% Repeatedly generate greedy solutions 
    stall_timer = tic;
    best_solution_set = {[] -Inf [] []}; % traj_set, demand_score, transfer_graph_digraph, solution,
    while toc(stall_timer) < overall_stall_time
        [solution, traj_set, score] = greedySolution(network, params, solution_stall_time);

        if best_solution_set{2} < score 
            stall_timer = tic;
            best_solution_set = {traj_set score solution};
            disp(strcat("Best score: ", string(best_solution_set{2})));
        end

        solution = best_solution_set{3};
        traj_set = best_solution_set{1};
    end
end

function [solution, traj_set] = geneticGlobalSearch(network, params)
    %% Run genetic algorithm globally over collision and objective
    nvars = params.n_trains * 2 * params.n_v_target_vars + params.n_switch_vars;
    % GA uses a nvars^2 matrix for mutation costing a lot of memory
    % GA wants normalized parameters
    obj_fun = @(solution) -combinedObjective(network, params, reshape(solution, params.n_trains, 2 * params.n_v_target_vars + params.n_switch_vars) + 0.5);
    options = optimoptions('ga', ...
        'Display','diagnose', ...
        'UseParallel', true, ...
        'MaxStallGenerations', 25, ...
        'PopulationSize', 15, ...
        'InitialPopulationRange', [-0.5; 0.5] ...
        );
    [X, fval, exitflag, output, population, scores] = ga(obj_fun, nvars, [], [], [], [], -0.5 * ones(nvars, 1), 0.5 * ones(nvars, 1), [], options);
    % Save result
    solution = reshape(X, params.n_trains, 2 * params.n_v_target_vars + params.n_switch_vars) + 0.5;
    [traj_set, ~] = constructTrajectorySet(network, params, solution);
end

function [solution, traj_set] = particleSwarmSearch(network, params)
    %% Run particle swarm optimization globally over collision and objective
    nvars = params.n_trains * 2 * params.n_v_target_vars + params.n_switch_vars;
    % GA uses a nvars^2 matrix for mutation costing a lot of memory
    % GA wants normalized parameters
    obj_fun = @(solution) -combinedObjective(network, params, reshape(solution, params.n_trains, 2 * params.n_v_target_vars + params.n_switch_vars) + 0.5);
    options = optimoptions('particleswarm', ...
        'Display','iter', ...
        'UseParallel', true, ...
        'MaxStallIterations', 25, ...
        'SwarmSize', 100, ...
        'InitialSwarmSpan', 1 ...
        );
    [X, fval, exitflag, output, scores] = particleswarm(obj_fun, nvars, -0.5 * ones(nvars, 1), 0.5 * ones(nvars, 1), options);
    % Save result
    solution = reshape(X, params.n_trains, 2 * params.n_v_target_vars + params.n_switch_vars) + 0.5;
    [traj_set, ~] = constructTrajectorySet(network, params, solution); 
end

function [solution, traj_set] = refineSolution(network, params, solution, traj_set)
    obj_fun = @(solution) -combinedObjective(network, params, reshape(solution, params.n_trains, 2 * params.n_v_target_vars + params.n_switch_vars));
    nvars = params.n_trains * 2 * params.n_v_target_vars + params.n_switch_vars;
    options = optimoptions('patternsearch', ...
        'Display','iter', ...
        'UseParallel', true, ...
        'InitialMeshSize', 0.1, ...
        'MeshTolerance', 0.02 ...
        );
    X = patternsearch(obj_fun, solution, [], [], [], [], zeros(nvars, 1), ones(nvars, 1), [], options);
    solution = reshape(X, params.n_trains, 2 * params.n_v_target_vars + params.n_switch_vars);
    [traj_set, ~, ~] = constructTrajectorySet(network, params, solution);
end