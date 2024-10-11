%% Generate Network
%network.adjacency_matrix = [ 0 500 200 0 500;
%                            0 0 600 2000 1000;
%                            0 0 0 400 0;
%                            0 0 0 0 0;
%                            0 0 0 0 0];
adj = random_planar_graph(3000);
network.adjacency_matrix = triu((adj + adj') * 1000, 1);
connection_indexes = find(network.adjacency_matrix!=0);
network.adjacency_matrix(connection_indexes) = network.adjacency_matrix(connection_indexes) .* randi([1,100], size(connection_indexes));
clear adj connection_indexes;
[network.edge_rows, network.edge_cols, network.edge_values] = find(network.adjacency_matrix);
network.adjacent_edge_list = {};
for node = 1:length(network.adjacency_matrix)(1)
    network.adjacent_edge_list{node} = find((network.edge_rows == node) | (network.edge_cols == node));
end

%% All shortest path pairs
distances = network.adjacency_matrix + network.adjacency_matrix';
distances(distances==0) = Inf;
distances(logical(eye(size(distances)))) = 0;
network.all_shortest_paths = FastFloyd(distances);
clear distances;

%% Parameters
params.n_timesteps = 8640; % 10s timesteps for one whole day
params.n_trains = 500;
params.min_separation = 100; % m
params.max_speed = 1.11; % m/10s 200km/h
params.max_accel = 46.27; % m/(10s)² 0-100kmh in 1m

% initial_speed dimensions (n_trains, 1) values (velocity 0-1)
params.initial_positions = randi([1,length(network.edge_values)], params.n_trains, 3);
params.initial_positions(:,2) = rand(params.n_trains,1);
params.initial_positions(:,3) = randi([0,1],params.n_trains,1) * 2 - 1;
params.initial_speeds = (rand(params.n_trains,1) * 2 - 1) * params.max_speed;

function traj_set = constructTrajectorySet(network, solution, initial_positions, initial_speeds, max_accel, max_speed)
    %% Constructs a set of train trajectories on the graph
    % solution dimensions (n_trains, timestep * 2)
    % solution values (acceleration 0-1, direction 0-1)
    % initial position dimensions (n_trains, 3)
    % initial position values (edge 0-n, position on edge 0-1)
    % initial speed dimensions (n_trains)
    % initial speed values (speed: -max_speed-max_speed)
    % trajectory dimensions (n_trains, 3, timestep)
    % trajectory values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)
    n_trains = size(solution)(1);
    for i_train = 1:n_trains
        traj_set(i_train, :, :) = constructTrajectory(network, solution(i_train,:), initial_positions(i_train, :), initial_speeds(i_train), max_accel, max_speed);
    end
end

function traj = constructTrajectory(network, solution, initial_position, initial_speed, max_accel, max_speed)
    %% Constructs a single train trajectory on the graph
    % solution dimensions (timestep * 2)
    % solution values (acceleration 0-1, direction 0-1)
    % initial position dimensions (3)
    % initial position values (edge 0-n, position on edge 0-1)
    % initial speed dimensions (1)
    % initial speed values (speed: -max_speed-max_speed)
    % trajectory dimensions (3, timestep)
    % trajectory values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)
    assert(~any(solution>1));
    assert(~any(solution<0));

    n_timesteps = size(solution)(2) / 2;

    %% Calculate speed curves
    % Speed is relative to train orientation
    acceleration = solution(1:n_timesteps);
    speeds = initial_speed + cumtrapz(((acceleration * 2) - 1) * max_accel);
    # Do not accelerate over maximum
    speeds(speeds>max_speed) = max_speed;
    speeds(speeds<-max_speed) = -max_speed;

    %% Calculate position curves
    position = cumtrapz(speeds);

    traj(1, :) = initial_position(1);
    traj(2, :) = initial_position(2);
    traj(3, :) = initial_position(3);

    pivot_timestep = 1;
    while pivot_timestep < n_timesteps
        current_edge_length = network.edge_values(traj(1, pivot_timestep));
        remaining_backward_length = traj(2, pivot_timestep) * current_edge_length;
        remaining_forward_length = current_edge_length - remaining_backward_length;
        base_position = position(pivot_timestep);

        % Find next edge change
        edge_trajectory = traj(3, pivot_timestep) * (position(pivot_timestep:end) - base_position);
        next_pivot_timestep = (pivot_timestep - 1) + find((edge_trajectory > remaining_forward_length | edge_trajectory < -remaining_backward_length), 1);
        forward_exit = any((edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) > remaining_forward_length));
        if isempty(next_pivot_timestep)
            next_pivot_timestep = n_timesteps;
        end
        assert(next_pivot_timestep <= n_timesteps);

        % Set trajectory on the current edge
        traj(1, pivot_timestep:next_pivot_timestep-1) = traj(1, pivot_timestep);
        traj(2, pivot_timestep:next_pivot_timestep-1) = traj(2, pivot_timestep) + edge_trajectory(1:next_pivot_timestep - (pivot_timestep-1) - 1) / current_edge_length;
        traj(3, pivot_timestep:next_pivot_timestep-1) = traj(3, pivot_timestep);

        %% Leave current edge
        if forward_exit
            traversed_node = network.edge_cols(traj(1, pivot_timestep));
            node_entrance_direction = traj(3, pivot_timestep);
            extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) - remaining_forward_length;
        else
            traversed_node = network.edge_rows(traj(1, pivot_timestep));
            node_entrance_direction = -traj(3, pivot_timestep);
            extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) + remaining_backward_length;
        end

        % Decide on next edge
        viable_next_edges = network.adjacent_edge_list{traversed_node};
        viable_next_edges = viable_next_edges(viable_next_edges!=traj(1, pivot_timestep));

        if isempty(viable_next_edges)
            % Stay stationary until trajectory comes back around
            # Use the old base position and forw/backw lengths
            deadend_edge_trajectory = traj(3, pivot_timestep) * (position(next_pivot_timestep:end) - base_position);
            if forward_exit
                deadend_next_pivot_timestep = (next_pivot_timestep - 1) + find(deadend_edge_trajectory < remaining_forward_length, 1);
            else
                deadend_next_pivot_timestep = (next_pivot_timestep - 1) + find(deadend_edge_trajectory > -remaining_backward_length, 1);
            end
            if isempty(deadend_next_pivot_timestep)
                deadend_next_pivot_timestep = n_timesteps;
            end
            assert(deadend_next_pivot_timestep <= n_timesteps);

            % Set position and train direction as stationary
            traj(1, next_pivot_timestep:deadend_next_pivot_timestep) = traj(1, pivot_timestep);
            traj(2, next_pivot_timestep:deadend_next_pivot_timestep) = double(forward_exit);
            traj(3, next_pivot_timestep:deadend_next_pivot_timestep) = traj(3, pivot_timestep);

            pivot_timestep = deadend_next_pivot_timestep;
        else
            next_edge_selection = 1 + round(solution(n_timesteps + pivot_timestep) * (length(viable_next_edges) - 1));
            next_edge = viable_next_edges(next_edge_selection);

            % Set position and train direction on new edge
            edge_entrance_direction = (network.edge_rows(next_edge) == traversed_node);
            traj(1, next_pivot_timestep) = next_edge; # edge number
            traj(2, next_pivot_timestep) = not(edge_entrance_direction) + (edge_entrance_direction*2 - 1) * abs(extra_movement / network.edge_values(next_edge)); # position on edge
            traj(3, next_pivot_timestep) = (edge_entrance_direction*2 - 1) * node_entrance_direction; # orientation on edge

            pivot_timestep = next_pivot_timestep;
        end
    end

    traj(:, n_timesteps) = traj(:, n_timesteps-1);
end

function score = objectiveFunction(network, traj_set, min_separation, max_speed)
    %% Evaluates a set of train trajectories
    % trajectory_set dimensions (n_trains, 3, timestep)
    % trajectory_set values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)
    penalty = 0;

    %% Separation Penalties
    % For each train pair update the minimum time to collision then skip that time and check again
    n_trains = size(traj_set, 1);
    n_timesteps = size(traj_set, 3);
    parfor i_train = 1:n_trains
        for j_train = i_train+1:n_trains
            timestep = 1;
            while timestep < n_timesteps
                distance = trainDistance(network, traj_set, i_train, j_train, timestep);

                if distance > min_separation
                    guaranteed_safe_time = int32(floor((distance - min_separation) / (2 * max_speed))) + 1;
                else
                    % Exponential penalty for closeness beyond the minimum separation
                    penalty = penalty + min(10000, (min_separation/distance - 1));
                    guaranteed_safe_time = 1;
                end

                timestep = timestep + guaranteed_safe_time;
            end
        end
    end

    %% Objective evaluation
    score = -penalty;
end

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
        %% Find shortest path with precomputed distance matrix
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

%% Solution Generation 

function solution = randomSolution(params)
    solution = rand(params.n_trains, params.n_timesteps * 2);
end

function solution = greedySolution(network, params)
    % Find valid routes sequentially
    solution = rand(1, params.n_timesteps * 2);
    traj_set(1,:,:) = constructTrajectory(network, solution, params.initial_positions(1,:), params.initial_speeds(1), params.max_accel, params.max_speed);
    i_train = 2;
    while i_train <= params.n_trains
        clc; disp(["Greedy Placement: ", mat2str(round(i_train/params.n_trains * 100)), "%"]);
        abort = 1;
        score = -Inf;
        while score < 0
            new_train_solution = rand(1, params.n_timesteps * 2);
            new_traj = constructTrajectory(network, new_train_solution, params.initial_positions(i_train,:), params.initial_speeds(i_train), params.max_accel, params.max_speed);
            new_traj_set = cat(1, traj_set, reshape(new_traj, 1, size(new_traj,1), size(new_traj,2)));
            new_solution = cat(1, solution, new_train_solution);
            score = objectiveFunction(network, new_traj_set, params.min_separation, params.max_speed);

            abort++;
            if abort > 1000
                warning("Failed to find valid placement. Trying from scratch..");
                new_solution = rand(1, params.n_timesteps * 2);
                new_traj_set = [];
                new_traj_set(1,:,:) = constructTrajectory(network, new_solution, params.initial_positions(1,:), params.initial_speeds(1), params.max_accel, params.max_speed);
                i_train = 1;
                break;
            end
        end
        traj_set = new_traj_set;
        solution = new_solution;
        i_train++;
    end
end

%% Search Methods

function randomSearch(network, params)
    csvwrite("network.csv", network.adjacency_matrix);
    score = -Inf;
    best_traj_set = [];
    while score < 0
        tic
        if greedy
            sol = greedySolution(network, params);
        else
        sol = randomSolution(params);
        end
        traj_set = constructTrajectorySet(network, sol, params.initial_positions, params.initial_speeds, params.max_accel, params.max_speed);
        new_score = objectiveFunction(network, traj_set, params.min_separation, params.max_speed);
        if score < new_score
            best_traj_set = traj_set;
            score = new_score
            csvwrite("trajectories_edges.csv", squeeze(traj_set(:,1,:)));
            csvwrite("trajectories_positions.csv", squeeze(traj_set(:,2,:)));
        end
    end
end

function localSearch(network, params)
    csvwrite("network.csv", network.adjacency_matrix);
    global_score = -Inf;

    for restart = 1:100
        improvement = Inf; 
        abort = 0;
        solution = randomSolution(params);
        score = objectiveFunction(network, constructTrajectorySet(network, solution, params.initial_positions, params.initial_speeds, params.max_accel, params.max_speed), params.min_separation, params.max_speed);

        while improvement > 100 && abort < 100
            preturbed_solution = solution + (rand(params.n_trains, params.n_timesteps * 2) - 0.5);
            preturbed_solution(preturbed_solution > 1) = 1;
            preturbed_solution(preturbed_solution < 0) = 0;
            traj_set = constructTrajectorySet(network, preturbed_solution, params.initial_positions, params.initial_speeds, params.max_accel, params.max_speed);
            new_score = objectiveFunction(network, traj_set, params.min_separation, params.max_speed);

            if new_score > score
                if new_score > global_score
                    csvwrite("trajectories_edges.csv", squeeze(traj_set(:,1,:)));
                    csvwrite("trajectories_positions.csv", squeeze(traj_set(:,2,:)));
                    global_score = new_score
                end
                solution = preturbed_solution;
                abort = 0;
                improvement = new_score - score;
            end
            abort++;
        end
    end
end

randomSearch(network, params, 1);