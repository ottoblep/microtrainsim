%% Generate Network
%params.adjacency_matrix = [ 0 500 200 0;
%                  0 0 600 0;
%                  0 0 0 400;
%                  0 0 0 0];
params.adjacency_matrix = triu(randomConnectedGraph(10,15) * 1000, 1);
[params.edge_rows, params.edge_cols, params.edge_values] = find(params.adjacency_matrix);
params.adjacent_edge_list = {};
for node = 1:length(params.adjacency_matrix)(1)
    params.adjacent_edge_list{node} = find((params.edge_rows == node) | (params.edge_cols == node));
end

%% All shortest path pairs
params.distances = params.adjacency_matrix + params.adjacency_matrix';
params.distances(params.distances==0) = Inf;
params.distances(logical(eye(size(params.distances)))) = 0;
params.all_shortest_paths = FastFloyd(params.distances);

%% Parameters
params.n_timesteps = 7640; % 10s timesteps for one whole day
params.n_trains = 10;
params.min_separation = 100; % m
params.max_speed = 1.11; % m/10s 200km/h
params.accel = 46.27; % m/(10s)² 0-100kmh in 1m

params.initial_pos = randi([1,length(params.edge_values)], params.n_trains, 3);
params.initial_pos(:,2) = rand(params.n_trains,1);
params.initial_pos(:,3) = randi([0,1],params.n_trains,1) * 2 - 1;
params.initial_speed = (rand(params.n_trains,1) * 2 - 1) * params.max_speed;

function sol = randomSolution(params)
    % Solution Array dimensions (train, timestep)
    % Solution Array values (acceleration -1,0,+1, direction 0-1)
    sol = rand(params.n_trains, params.n_timesteps, 2);
    sol(:,:,1) = randi([-1,1], params.n_trains, params.n_timesteps);
end

function traj = constructTrajectory(params, solution)
    % solution dimensions (train, timestep) values (acceleration -1,0,+1, direction 0-1)
    % initial_pos dimensions (n_trains, 2) values (edge 1-n, position on edge 0-1
    % initial_speed dimensions (n_trains, 1) values (velocity 0-1)
    % trajectory dimensions (train, timestep) values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)

    %% Calculate speed curves
    % Speed is relative to train orientation
    parfor i_train = 1:params.n_trains
        speeds(i_train, :) = params.initial_speed(i_train) + cumtrapz(solution(i_train, :, 1)) * params.accel;
        # Do not accelerate over maximum
        speeds(i_train, speeds(i_train,:)>params.max_speed) = params.max_speed;
        speeds(i_train,speeds(i_train,:)<-params.max_speed) = -params.max_speed;
    end

    %% Calculate position curves
    parfor i_train = 1:params.n_trains
        position(i_train, :) = cumtrapz(speeds(i_train, :));
    end

    traj(:, 1, :) = params.initial_pos;

    parfor i_train = 1:params.n_trains
        pivot_timestep = 1;
        while pivot_timestep < params.n_timesteps
            current_edge_length = params.edge_values(traj(i_train, pivot_timestep, 1));
            remaining_backward_length = traj(i_train, pivot_timestep, 2) * current_edge_length;
            remaining_forward_length = current_edge_length - remaining_backward_length;
            base_position = position(i_train, pivot_timestep);

            % Find next edge change
            edge_trajectory = traj(i_train, pivot_timestep, 3) * (position(i_train,pivot_timestep:end) - base_position);
            next_pivot_timestep = (pivot_timestep - 1) + find((edge_trajectory > remaining_forward_length | edge_trajectory < -remaining_backward_length), 1);
            forward_exit = any((edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) > remaining_forward_length));
            if isempty(next_pivot_timestep)
                next_pivot_timestep = params.n_timesteps;
            end
            assert(next_pivot_timestep <= params.n_timesteps);

            % Set trajectory on the current edge
            traj(i_train, pivot_timestep:next_pivot_timestep-1, 1) = traj(i_train, pivot_timestep, 1);
            traj(i_train, pivot_timestep:next_pivot_timestep-1, 2) = traj(i_train, pivot_timestep, 2) + edge_trajectory(1:next_pivot_timestep - (pivot_timestep-1) - 1) / current_edge_length;
            traj(i_train, pivot_timestep:next_pivot_timestep-1, 3) = traj(i_train, pivot_timestep, 3);

            %% Leave current edge
            if forward_exit
                traversed_node = params.edge_cols(traj(i_train, pivot_timestep, 1));
                node_entrance_direction = traj(i_train, pivot_timestep, 3);
                extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) - remaining_forward_length;
            else
                traversed_node = params.edge_rows(traj(i_train, pivot_timestep, 1));
                node_entrance_direction = -traj(i_train, pivot_timestep, 3);
                extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) + remaining_backward_length;
            end

            % Decide on next edge
            viable_next_edges = params.adjacent_edge_list{traversed_node};
            viable_next_edges = viable_next_edges(viable_next_edges!=traj(i_train, pivot_timestep, 1));

            if isempty(viable_next_edges)
                % Stay stationary until trajectory comes back around
                # Use the old base position and forw/backw lengths
                deadend_edge_trajectory = traj(i_train, pivot_timestep, 3) * (position(i_train,next_pivot_timestep:end) - base_position);
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
                traj(i_train, next_pivot_timestep:deadend_next_pivot_timestep, 1) = traj(i_train, pivot_timestep, 1);
                traj(i_train, next_pivot_timestep:deadend_next_pivot_timestep, 2) = double(forward_exit);
                traj(i_train, next_pivot_timestep:deadend_next_pivot_timestep, 3) = traj(i_train, pivot_timestep, 3);

                pivot_timestep = deadend_next_pivot_timestep;
            else
                next_edge = viable_next_edges(int32(floor(solution(i_train, pivot_timestep, 2) * length(viable_next_edges)) + 1));

                % Set position and train direction on new edge
                edge_entrance_direction = (params.edge_rows(next_edge) == traversed_node);
                traj(i_train, next_pivot_timestep, 1) = next_edge; # edge number
                traj(i_train,next_pivot_timestep, 2) = not(edge_entrance_direction) + (edge_entrance_direction*2 - 1) * abs(extra_movement / params.edge_values(next_edge)); # position on edge
                traj(i_train, next_pivot_timestep, 3) = (edge_entrance_direction*2 - 1) * node_entrance_direction; # orientation on edge

                pivot_timestep = next_pivot_timestep;
            end
        end

        traj(i_train, params.n_timesteps, :) = traj(i_train, params.n_timesteps-1, :);
    end
end

function score = objectiveFunction(params, traj)
    penalty = 0;

    %% Separation Penalties
    % For each train pair update the minimum time to collision then skip that time and check again
    parfor i_train = 1:params.n_trains
        timestep = 1;
        for j_train = i_train+1:params.n_trains
            while timestep < params.n_timesteps
                edge_i = int32(traj(i_train, timestep, 1));
                edge_j = int32(traj(j_train, timestep, 1));
                i_edge_length = params.edge_values(edge_i);
                j_edge_length = params.edge_values(edge_j);

                if edge_i == edge_j
                    distance = abs(traj(i_train, timestep, 2) - traj(j_train, timestep, 2)) * i_edge_length;
                else
                    %% Find shortest path with precomputed distance matrix
                    i_node_backward = params.edge_rows(edge_i);
                    i_node_forward = params.edge_cols(edge_i);
                    j_node_backward = params.edge_rows(edge_j);
                    j_node_forward = params.edge_cols(edge_j);

                    i_remaining_backward_length = traj(i_train, timestep, 2) * i_edge_length;
                    i_remaining_forward_length = i_edge_length - i_remaining_backward_length;
                    j_remaining_backward_length = traj(j_train, timestep, 2) * j_edge_length;
                    j_remaining_forward_length = j_edge_length - j_remaining_backward_length;

                    dist1 = i_remaining_backward_length + params.all_shortest_paths(i_node_backward,j_node_backward) + j_remaining_backward_length;
                    dist2 = i_remaining_forward_length + params.all_shortest_paths(i_node_forward,j_node_backward) + j_remaining_backward_length;
                    dist3 = i_remaining_backward_length + params.all_shortest_paths(i_node_backward,j_node_forward) + j_remaining_forward_length;
                    dist4 = i_remaining_forward_length + params.all_shortest_paths(i_node_forward,j_node_forward) + j_remaining_forward_length;

                    distance = min([dist1, dist2, dist3, dist4]);
                end

                if distance > params.min_separation
                    guaranteed_safe_time = int32(floor((distance - params.min_separation) / (2 * params.max_speed))) + 1;
                else
                    % Exponential penalty for closeness beyond the minimum separation
                    penalty = penalty + (params.min_separation/distance - 1);
                    guaranteed_safe_time = 1;
                end

                timestep = timestep + guaranteed_safe_time;
            end
        end
    end

    %% Objective evaluation
    score = -penalty;
end

tic
sol = randomSolution(params);
traj = constructTrajectory(params, sol);
toc
tic
score = objectiveFunction(params, traj)
toc
csvwrite("network.csv", params.adjacency_matrix);
csvwrite("trajectories_edges.csv", traj(:,:,1));
csvwrite("trajectories_positions.csv", traj(:,:,2));
