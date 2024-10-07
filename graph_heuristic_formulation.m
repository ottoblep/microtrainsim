%% Generate Network
%params.adjacency_matrix = [ 0 500 200 0 500;
%                            0 0 600 2000 1000;
%                            0 0 0 400 0;
%                            0 0 0 0 0;
%                            0 0 0 0 0];
adj = randomGraph(10,20);
params.adjacency_matrix = triu((adj + adj') * 1000, 1);
clear adj;
[params.edge_rows, params.edge_cols, params.edge_values] = find(params.adjacency_matrix);
params.adjacent_edge_list = {};
for node = 1:length(params.adjacency_matrix)(1)
    params.adjacent_edge_list{node} = find((params.edge_rows == node) | (params.edge_cols == node));
end

%% All shortest path pairs
distances = params.adjacency_matrix + params.adjacency_matrix';
distances(distances==0) = Inf;
distances(logical(eye(size(distances)))) = 0;
params.all_shortest_paths = FastFloyd(distances);
clear distances;

%% Parameters
params.n_timesteps = 8640; % 10s timesteps for one whole day
params.n_trains = 10;
params.min_separation = 100; % m
params.max_speed = 1.11; % m/10s 200km/h
params.max_accel = 46.27; % m/(10s)² 0-100kmh in 1m

% initial_speed dimensions (n_trains, 1) values (velocity 0-1)
params.initial_positions = randi([1,length(params.edge_values)], params.n_trains, 3);
params.initial_positions(:,2) = rand(params.n_trains,1);
params.initial_positions(:,3) = randi([0,1],params.n_trains,1) * 2 - 1;
params.initial_speeds = (rand(params.n_trains,1) * 2 - 1) * params.max_speed;

function traj_set = constructTrajectorySet(params, solution, initial_positions, initial_speeds)
    %% Constructs a set of train trajectories on the graph
    % solution dimensions (timestep * 2)
    % solution values (acceleration 0-1, direction 0-1)
    % initial position dimensions (n_trains, 3)
    % initial position values (edge 0-n, position on edge 0-1)
    % initial speed dimensions (n_trains)
    % initial speed values (speed: -max_speed-max_speed)
    % trajectory dimensions (n_trains, 3, timestep)
    % trajectory values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)
    for i_train = 1:params.n_trains
        accelerations = solution((i_train - 1) * params.n_timesteps + 1:(i_train + 1) * params.n_timesteps);
        directions = solution((params.n_trains + i_train - 1) * params.n_timesteps + 1:(params.n_trains + i_train) * params.n_timesteps);
        single_solution = cat(2, accelerations, directions);
        traj_set(i_train, :, :) = constructTrajectory(params, single_solution, initial_positions(i_train, :), initial_speeds(i_train));
    end
end

function traj = constructTrajectory(params, solution, initial_position, initial_speed)
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

    %% Calculate speed curves
    % Speed is relative to train orientation
    acceleration = solution(1:params.n_timesteps);
    speeds = initial_speed + cumtrapz(((acceleration * 2) - 1) * params.max_accel);
    # Do not accelerate over maximum
    speeds(speeds>params.max_speed) = params.max_speed;
    speeds(speeds<-params.max_speed) = -params.max_speed;

    %% Calculate position curves
    position = cumtrapz(speeds);

    traj(1, :) = initial_position(1);
    traj(2, :) = initial_position(2);
    traj(3, :) = initial_position(3);

    pivot_timestep = 1;
    while pivot_timestep < params.n_timesteps
        current_edge_length = params.edge_values(traj(1, pivot_timestep));
        remaining_backward_length = traj(2, pivot_timestep) * current_edge_length;
        remaining_forward_length = current_edge_length - remaining_backward_length;
        base_position = position(pivot_timestep);

        % Find next edge change
        edge_trajectory = traj(3, pivot_timestep) * (position(pivot_timestep:end) - base_position);
        next_pivot_timestep = (pivot_timestep - 1) + find((edge_trajectory > remaining_forward_length | edge_trajectory < -remaining_backward_length), 1);
        forward_exit = any((edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) > remaining_forward_length));
        if isempty(next_pivot_timestep)
            next_pivot_timestep = params.n_timesteps;
        end
        assert(next_pivot_timestep <= params.n_timesteps);

        % Set trajectory on the current edge
        traj(1, pivot_timestep:next_pivot_timestep-1) = traj(1, pivot_timestep);
        traj(2, pivot_timestep:next_pivot_timestep-1) = traj(2, pivot_timestep) + edge_trajectory(1:next_pivot_timestep - (pivot_timestep-1) - 1) / current_edge_length;
        traj(3, pivot_timestep:next_pivot_timestep-1) = traj(3, pivot_timestep);

        %% Leave current edge
        if forward_exit
            traversed_node = params.edge_cols(traj(1, pivot_timestep));
            node_entrance_direction = traj(3, pivot_timestep);
            extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) - remaining_forward_length;
        else
            traversed_node = params.edge_rows(traj(1, pivot_timestep));
            node_entrance_direction = -traj(3, pivot_timestep);
            extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) + remaining_backward_length;
        end

        % Decide on next edge
        viable_next_edges = params.adjacent_edge_list{traversed_node};
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
                deadend_next_pivot_timestep = params.n_timesteps;
            end
            assert(deadend_next_pivot_timestep <= params.n_timesteps);

            % Set position and train direction as stationary
            traj(1, next_pivot_timestep:deadend_next_pivot_timestep) = traj(1, pivot_timestep);
            traj(2, next_pivot_timestep:deadend_next_pivot_timestep) = double(forward_exit);
            traj(3, next_pivot_timestep:deadend_next_pivot_timestep) = traj(3, pivot_timestep);

            pivot_timestep = deadend_next_pivot_timestep;
        else
            next_edge_selection = 1 + round(solution(params.n_timesteps + pivot_timestep) * (length(viable_next_edges) - 1));
            next_edge = viable_next_edges(next_edge_selection);

            % Set position and train direction on new edge
            edge_entrance_direction = (params.edge_rows(next_edge) == traversed_node);
            traj(1, next_pivot_timestep) = next_edge; # edge number
            traj(2, next_pivot_timestep) = not(edge_entrance_direction) + (edge_entrance_direction*2 - 1) * abs(extra_movement / params.edge_values(next_edge)); # position on edge
            traj(3, next_pivot_timestep) = (edge_entrance_direction*2 - 1) * node_entrance_direction; # orientation on edge

            pivot_timestep = next_pivot_timestep;
        end
    end

    traj(:, params.n_timesteps) = traj(:, params.n_timesteps-1);
end

function score = objectiveFunction(params, traj_set)
    %% Evaluates a set of train trajectories
    % trajectory_set dimensions (n_trains, 3, timestep)
    % trajectory_set values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)
    penalty = 0;

    %% Separation Penalties
    % For each train pair update the minimum time to collision then skip that time and check again
    parfor i_train = 1:params.n_trains
        for j_train = i_train+1:params.n_trains
            timestep = 1;
            while timestep < params.n_timesteps
                distance = trainDistance(params, traj_set, i_train, j_train, timestep);

                if distance > params.min_separation
                    guaranteed_safe_time = int32(floor((distance - params.min_separation) / (2 * params.max_speed))) + 1;
                else
                    % Exponential penalty for closeness beyond the minimum separation
                    penalty = penalty + min(10000, (params.min_separation/distance - 1));
                    guaranteed_safe_time = 1;
                end

                timestep = timestep + guaranteed_safe_time;
            end
        end
    end

    %% Objective evaluation
    score = -penalty;
end

function distance = trainDistance(params, traj_set, i_train, j_train, timestep)
    %% Calculates distance of two trains given a set of trajectories
    % trajectory_set dimensions (n_trains, 3, timestep)
    % trajectory_set values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)

    edge_i = int32(traj_set(i_train, 1, timestep));
    edge_j = int32(traj_set(j_train, 1, timestep));

    i_edge_length = params.edge_values(edge_i);
    j_edge_length = params.edge_values(edge_j);

    if edge_i == edge_j
        distance = abs(traj_set(i_train, 2, timestep) - traj_set(j_train, 2, timestep)) * i_edge_length;
    else
        %% Find shortest path with precomputed distance matrix
        i_node_backward = params.edge_rows(edge_i);
        i_node_forward = params.edge_cols(edge_i);
        j_node_backward = params.edge_rows(edge_j);
        j_node_forward = params.edge_cols(edge_j);

        i_remaining_backward_length = traj_set(i_train, 2, timestep) * i_edge_length;
        i_remaining_forward_length = i_edge_length - i_remaining_backward_length;
        j_remaining_backward_length = traj_set(j_train, 2, timestep) * j_edge_length;
        j_remaining_forward_length = j_edge_length - j_remaining_backward_length;

        dist1 = i_remaining_backward_length + params.all_shortest_paths(i_node_backward,j_node_backward) + j_remaining_backward_length;
        dist2 = i_remaining_forward_length + params.all_shortest_paths(i_node_forward,j_node_backward) + j_remaining_backward_length;
        dist3 = i_remaining_backward_length + params.all_shortest_paths(i_node_backward,j_node_forward) + j_remaining_forward_length;
        dist4 = i_remaining_forward_length + params.all_shortest_paths(i_node_forward,j_node_forward) + j_remaining_forward_length;

        distance = min([dist1, dist2, dist3, dist4]);
    end
end

%% Search algorithms

function sol = greedySolution(params)
end

function sol = randomSolution(params)
    sol = rand(1, params.n_trains * params.n_timesteps * 2);
end

function randomSearch(params)
    csvwrite("network.csv", params.adjacency_matrix);
    score = -Inf;
    best_traj_set = [];
    while score < 0
        tic
        sol = randomSolution(params);
        traj_set = constructTrajectorySet(params, sol, params.initial_positions, params.initial_speeds);
        new_score = objectiveFunction(params, traj_set);
        if score < new_score
            best_traj_set = traj_set;
            score = new_score
            csvwrite("trajectories_edges.csv", squeeze(traj_set(:,1,:)));
            csvwrite("trajectories_positions.csv", squeeze(traj_set(:,2,:)));
        end
    end
end

function localSearch(params)
    csvwrite("network.csv", params.adjacency_matrix);
    global_score = -Inf;

    for restart = 1:100
        improvement = Inf; 
        abort = 0;
        solution = randomSolution(params);
        score = objectiveFunction(params,constructTrajectorySet(params, solution, params.initial_positions, params.initial_speeds));

        while improvement > 100 && abort < 100
            preturbed_solution = solution + (rand(1, params.n_trains * params.n_timesteps * 2) - 0.5);
            preturbed_solution(preturbed_solution > 1) = 1;
            preturbed_solution(preturbed_solution < 0) = 0;
            traj_set = constructTrajectorySet(params, preturbed_solution, params.initial_positions, params.initial_speeds);
            new_score = objectiveFunction(params, traj_set);

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

localSearch(params);
