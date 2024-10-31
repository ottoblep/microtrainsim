function greedyHeuristicFormulation()
    [network, params] = generateEnvironment();

    csvwrite("network.csv", network.adjacency_matrix);
    csvwrite("demands.csv", params.demand_matrix);

    %[solution, traj_set] = greedySearch(network, params, true, 2, 2);
    %disp("---");
    %[solution, traj_set] = greedySearch(network, params, false, 10, 0.01);
    %[solution, traj_set] = geneticGlobalSearch(network, params);
    [solution, traj_set] = particleSwarmSearch(network, params);

    csvwrite("trajectories_edges.csv", squeeze(traj_set(:,1,:)));
    csvwrite("trajectories_positions.csv", squeeze(traj_set(:,2,:)));
end

function [network, params] = generateEnvironment()
    %% Network
    adj = randomPlanarGraph(5);
    network.adjacency_matrix = triu((adj + adj') * 1000, 1);
    network.adjacency_matrix = [ 0 0 0 0 0 0 1000 0 0 0 0 0 0 0 0;
                                 0 0 0 0 0 0 1000 1000 0 0 0 0 0 0 0;
                                 0 0 0 0 0 0 0 1000 0 0 0 0 0 0 0;
                                 0 0 0 0 0 0 1000 0 0 0 0 0 0 0 0;
                                 0 0 0 0 0 0 0 1000 0 0 0 0 0 0 0;
                                 0 0 0 0 0 0 1000 0 0 0 0 0 0 0 0; 
                                 0 0 0 0 0 0 0 1000 0 1000 1000 0 0 0 0;
                                 0 0 0 0 0 0 0 0 1000 0 1000 1000 0 0 0;
                                 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                 0 0 0 0 0 0 0 0 0 0 1000 0 0 0 0;
                                 0 0 0 0 0 0 0 0 0 0 0 1000 1000 1000 1000;
                                 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
                                 ];
    connection_indexes = find(network.adjacency_matrix);
    network.adjacency_matrix(connection_indexes) = network.adjacency_matrix(connection_indexes) .* (randi([1,3],size(connection_indexes)));
    [network.edge_rows, network.edge_cols, network.edge_values] = find(network.adjacency_matrix);

    network.adjacent_edge_list = {};
    for node = 1:size(network.adjacency_matrix,1)
        network.adjacent_edge_list{node} = find((network.edge_rows == node) | (network.edge_cols == node));
    end

    %% All shortest path pairs
    tmp_adj = network.adjacency_matrix;
    tmp_adj(tmp_adj==0) = Inf;
    network.all_shortest_paths = distances(graph(tmp_adj,'upper'), 'Method', 'positive');

    %% Parameters
    params.n_timesteps = 8640; % 10s timesteps for one whole day, must be divisible my interpolation factor
    params.interpolation_factor = 160; % Support points for acceleration and switch direction curves for every n timesteps (not linearly spaced)
    params.n_trains = 12;
    params.min_separation = 100; % m
    params.max_speed = 1.11; % m/10s 200km/h
    params.max_accel = 46.27; % m/(10s)² 0-100kmh in 1m
    params.max_changeover_time = 1440; % 4hrs
    params.train_capacity = 400; % 400 Passengers

    assert(size(network.adjacency_matrix,1) > params.n_trains); % Trains can start at a unique node
    assert(mod(params.n_timesteps, params.interpolation_factor) == 0);

    params.demand_matrix = randi(1000, size(network.adjacency_matrix,1));
    params.demand_matrix(logical(eye(size(params.demand_matrix)))) = 0;
    params.demand_matrix(randperm(numel(params.demand_matrix), round(numel(params.demand_matrix)/2))) = 0; % make demand matrix more sparse
    %params.demand_matrix = [ 0 20;
    %                         0 0];

    % Random initial position
    %params.initial_positions = randi([1,length(network.edge_values)], params.n_trains, 3);
    %params.initial_positions(:,2) = rand(params.n_trains,1);
    %params.initial_positions(:,3) = randi([0,1],params.n_trains,1) * 2 - 1;
    %params.initial_speeds = (rand(params.n_trains,1) * 2 - 1) * params.max_speed;
    %params.destinations = randi([1,size(network.adjacency_matrix)], params.n_trains, 1);

    % Unique node initial position
    start_nodes = [1 2 3 4 6 10 5 9 12 13 14 15];
    for i_train = 1:params.n_trains
        [params.initial_positions(i_train, 1), params.initial_positions(i_train, 2)] = nodeToEdgePos(network, start_nodes(i_train));
    end
    params.initial_positions(:, 3) = ones(params.n_trains, 1);
    params.initial_speeds = zeros(params.n_trains, 1);

    params.destinations = [13 14 15 5 9 12 4 6 10 1 2 3];
    assert(numel(params.destinations) == params.n_trains);
end

%% Solution Construction

function [traj_set, event_set] = constructTrajectorySet(network, solution, initial_positions, initial_speeds, max_accel, max_speed, interpolation_factor)
    %% Constructs a set of train trajectories on the graph
    % solution dimensions (n_trains, 4 * timestep / interpolation_factor)
    % initial position dimensions (n_trains, 3)
    % initial speed dimensions (n_trains)
    % traj_set dimensions (n_trains, 3, timestep)
    % event_set dimensions (3, n))
    % event_set values (presence_at_node, timestep, train)
    n_trains = size(solution,1);
    event_set = [];
    traj_set = zeros(n_trains, 3, size(solution,2) * interpolation_factor / 4);
    for i_train = 1:n_trains
        [traj_set(i_train, :, :), new_events] = constructTrajectory(network, solution(i_train,:), initial_positions(i_train, :), initial_speeds(i_train), max_accel, max_speed, interpolation_factor);
        new_events(3,:) = i_train;
        event_set = cat(2, event_set, new_events);
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
    %% Penalizes train's distance from destination at the end of the trajectory
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
    %% Evaluates satisfied demand for a given set of arrival/departure timesteps 
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
    [traj_set, event_set] = constructTrajectorySet(network, solution, params.initial_positions, params.initial_speeds, params.max_accel, params.max_speed, params.interpolation_factor);
    collision_score = collisionPenalties(network, traj_set, params.min_separation, params.max_speed);
    [demand_score, ~, ~] = demandSatisfaction(network, event_set, params.demand_matrix, params.max_changeover_time, params.train_capacity);
    destination_score = destinationPenalties(network, traj_set, params.destinations);
    combined_score = collision_score + demand_score + destination_score;
end

%% Solution Generation 

function solution = randomSolution(params)
    solution = rand(params.n_trains, 4 * params.n_timesteps / params.interpolation_factor);
end

function solution = greedyValidSolution(network, params, stall_time)
    %% Brute-forces valid routes sequentially
    solution = rand(1, 4 * params.n_timesteps / params.interpolation_factor);;
    traj_set(1,:,:) = constructTrajectory(network, solution, params.initial_positions(1,:), params.initial_speeds(1), params.max_accel, params.max_speed, params.interpolation_factor);
    i_train = 2;
    while i_train <= params.n_trains
        stall_timer = tic;
        collision_score = -Inf;
        while collision_score < 0
            new_train_solution = rand(1, 4 * params.n_timesteps / params.interpolation_factor);
            new_traj = constructTrajectory(network, new_train_solution, params.initial_positions(i_train,:), params.initial_speeds(i_train), params.max_accel, params.max_speed, params.interpolation_factor);
            new_traj_set = cat(1, traj_set, reshape(new_traj, 1, size(new_traj,1), size(new_traj,2)));
            new_solution = cat(1, solution, new_train_solution);
            collision_score = collisionPenalties(network, new_traj_set, params.min_separation, params.max_speed);

            if toc(stall_timer) > stall_time
                disp("Failed to construct a collision free solution. Restarting ...");
                new_solution = rand(1, 4 * params.n_timesteps / params.interpolation_factor);
                new_traj_set = [];
                new_traj_set(1,:,:) = constructTrajectory(network, new_train_solution, params.initial_positions(1,:), params.initial_speeds(1), params.max_accel, params.max_speed, params.interpolation_factor);
                i_train = 1;
                break;
            end
        end
        traj_set = new_traj_set;
        solution = new_solution;
        i_train = i_train + 1;
    end
end

function [solution, traj_set] = greedySolution(network, params, per_train_stall_time)
    %% Place each train greedily considering demand and collisions
    traj_set = [];
    event_set = [];
    solution = [];

    for i_train = 1:params.n_trains
        stall_timer = tic;
        round_best_traj_set = [];
        round_best_event_set = [];
        round_best_train_solution = [];
        round_best_score = -Inf;

        while toc(stall_timer) < per_train_stall_time
            % Generate new trajectory
            new_train_solution = rand(1, 4 * params.n_timesteps / params.interpolation_factor);
            [new_traj, new_events] = constructTrajectory(network, new_train_solution, params.initial_positions(i_train,:), params.initial_speeds(i_train), params.max_accel, params.max_speed, params.interpolation_factor);

            % Add new trajectory to old set
            new_traj_set = cat(1, traj_set, reshape(new_traj, 1, size(new_traj,1), size(new_traj,2)));
            new_events(3,:) = i_train;
            new_event_set = cat(2, event_set, new_events);

            % Test new trajectory set
            collision_score = collisionPenalties(network, new_traj_set, params.min_separation, params.max_speed);
            demand_score = demandSatisfaction(network, event_set, params.demand_matrix, params.max_changeover_time, params.train_capacity);
            destination_score = destinationPenalties(network, new_traj_set, params.destinations);

            if collision_score + demand_score + destination_score > round_best_score
                round_best_score = collision_score + demand_score + destination_score;
                stall_timer = tic;
                round_best_event_set = new_event_set;
                round_best_traj_set = new_traj_set;
                round_best_train_solution = new_train_solution;
            end
        end

        traj_set = round_best_traj_set;
        event_set = round_best_event_set;
        solution = cat(1, solution, round_best_train_solution);
    end
end

%% Search Methods

function [solution, traj_set] = greedySearch(network, params, valid_solutions_suffice, overall_stall_time, solution_stall_time)
    %% Generate greedy solutions then evaluate demand score
    stall_timer = tic;
    best_solution_set = {[] -Inf [] []}; % traj_set, demand_score, transfer_graph_digraph, solution,
    while toc(stall_timer) < overall_stall_time
        if valid_solutions_suffice
            solution = greedyValidSolution(network, params, solution_stall_time);
        else
            solution = greedySolution(network, params, solution_stall_time);
        end
        [traj_set, event_set] = constructTrajectorySet(network, solution, params.initial_positions, params.initial_speeds, params.max_accel, params.max_speed, params.interpolation_factor);
        collision_score = collisionPenalties(network, traj_set, params.min_separation, params.max_speed);
        [demand_score, ~, ~] = demandSatisfaction(network, event_set, params.demand_matrix, params.max_changeover_time, params.train_capacity);
        destination_score = destinationPenalties(network, traj_set, params.destinations);
        new_score = collision_score + demand_score + destination_score;

        if best_solution_set{2} < new_score 
            stall_timer = tic;
            best_solution_set = {traj_set new_score solution};
            disp(strcat("Best score: ", string(best_solution_set{2})));
        end

        solution = best_solution_set{3};
        traj_set = best_solution_set{1};
    end
end

function [solution, traj_set] = geneticGlobalSearch(network, params)
    %% Run genetic algorithm globally over collision and demand objective
    nvars = params.n_trains * 4 * params.n_timesteps / params.interpolation_factor;
    % GA uses a nvars^2 matrix for mutation costing a lot of memory
    % GA wants normalized parameters
    obj_fun = @(solution) -combinedObjective(network, params, reshape(solution, params.n_trains, 4 * params.n_timesteps / params.interpolation_factor) + 0.5);
    options = optimoptions('ga', ...
        'Display','diagnose', ...
        'UseParallel', true, ...
        'MaxStallGenerations', 25, ...
        'PopulationSize', 15, ...
        'InitialPopulationRange', [-0.5; 0.5] ...
        );
    [X, fval, exitflag, output, population, scores] = ga(obj_fun, nvars, [], [], [], [], -0.5 * ones(nvars, 1), 0.5 * ones(nvars, 1), [], options);
    % Save result
    solution = reshape(X, params.n_trains, 4 * params.n_timesteps / params.interpolation_factor) + 0.5;
    [traj_set, event_set] = constructTrajectorySet(network, solution, params.initial_positions, params.initial_speeds, params.max_accel, params.max_speed, params.interpolation_factor);
end

function [solution, traj_set] = particleSwarmSearch(network, params)
    %% Run particle swarm optimization globally over collision and demand objective
    nvars = params.n_trains * 4 * params.n_timesteps / params.interpolation_factor;
    % GA uses a nvars^2 matrix for mutation costing a lot of memory
    % GA wants normalized parameters
    obj_fun = @(solution) -combinedObjective(network, params, reshape(solution, params.n_trains, 4 * params.n_timesteps / params.interpolation_factor) + 0.5);
    options = optimoptions('particleswarm', ...
        'Display','iter', ...
        'UseParallel', true, ...
        'MaxStallIterations', 25, ...
        'SwarmSize', 100, ...
        'InitialSwarmSpan', 1 ...
        );
    [X, fval, exitflag, output, scores] = particleswarm(obj_fun, nvars, -0.5 * ones(nvars, 1), 0.5 * ones(nvars, 1), options);
    % Save result
    solution = reshape(X, params.n_trains, 4 * params.n_timesteps / params.interpolation_factor) + 0.5;
    [traj_set, event_set] = constructTrajectorySet(network, solution, params.initial_positions, params.initial_speeds, params.max_accel, params.max_speed, params.interpolation_factor);
end

function [solution, traj_set] = refineSolution(network, params, solution, traj_set)
    obj_fun = @(solution) -combinedObjective(network, params, reshape(solution, params.n_trains, 4 * params.n_timesteps / params.interpolation_factor));
    nvars = params.n_trains * 4 * params.n_timesteps / params.interpolation_factor;
    options = optimoptions('patternsearch', ...
        'Display','iter', ...
        'UseParallel', true, ...
        'InitialMeshSize', 0.1, ...
        'MeshTolerance', 0.02 ...
        );
    X = patternsearch(obj_fun, solution, [], [], [], [], zeros(nvars, 1), ones(nvars, 1), [], options);
    solution = reshape(X, params.n_trains, 4 * params.n_timesteps / params.interpolation_factor);
    [traj_set, ~] = constructTrajectorySet(network, solution, params.initial_positions, params.initial_speeds, params.max_accel, params.max_speed, params.interpolation_factor);
end

function [solution, traj_set] = repairHeuristic(network, params, solution, traj_set)
    %% Preturb a solution until collisions disappear
    n_trains = size(traj_set, 1);
    n_timesteps = size(traj_set, 3);
    flawless = false;
    
    while ~flawless
        flawless = true;
        for i_train = 1:n_trains
            for j_train = i_train+1:n_trains
                timestep = 1;
                while timestep < n_timesteps
                    distance = trainDistance(network, traj_set, i_train, j_train, timestep);

                    if distance > params.min_separation
                        guaranteed_safe_time = int32(floor((distance - params.min_separation) / (2 * params.max_speed))) + 1;
                    else
                        solution(i_train, :) = solution(i_train, :) + rand(1, 4 * params.n_timesteps / params.interpolation_factor) * 0.05 - 0.025;
                        solution(i_train, solution(i_train, :) > 1) = 1;
                        solution(i_train, solution(i_train, :) < 0) = 0;
                        solution(j_train, :) = solution(i_train, :) + rand(1, 4 * params.n_timesteps / params.interpolation_factor) * 0.05 - 0.025;
                        solution(j_train, solution(j_train, :) > 1) = 1;
                        solution(j_train, solution(j_train, :) < 0) = 0;

                        [traj_set(i_train, :, :), ~] = constructTrajectory(network, solution(i_train,:), params.initial_positions(i_train, :), params.initial_speeds(i_train), params.max_accel, params.max_speed, params.interpolation_factor);
                        [traj_set(j_train, :, :), ~] = constructTrajectory(network, solution(j_train,:), params.initial_positions(j_train, :), params.initial_speeds(j_train), params.max_accel, params.max_speed, params.interpolation_factor);

                        flawless = false;
                        timestep =  1;
                        continue;
                    end

                    timestep = timestep + guaranteed_safe_time;
                end
            end
        end
    end

    collision_score = collisionPenalties(network, traj_set, params.min_separation, params.max_speed)
end