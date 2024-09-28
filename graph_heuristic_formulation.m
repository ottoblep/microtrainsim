params.network = [ 0 50 10 0;
                 0 0 60 0;
                 0 0 0 40;
                 0 0 0 0];
[params.edge_rows, params.edge_cols, params.edge_values] = find(params.network);
params.distances = params.network;
params.distances(params.network==0) = Inf;
params.all_shortest_paths = FastFloyd(params.distances);

params.n_trains = 3;
params.initial_pos = [2, 0.8, 1; 3, 0.7, 1; 4, 0.01, -1;];
params.initial_speed = [0.2, 0.3, 0.4];
params.n_timesteps = 1000;
params.min_separation = 5;
params.max_speed = 1;
params.accel = 0.1;

function sol = randomSolution(params)
    % Solution Array dimensions (train, timestep)
    % Solution Array values (acceleration -1,0,+1, direction 0-1)
    sol = rand(params.n_trains, params.n_timesteps - 1, 2);
    sol(:,:,1) = randi([-1,1], params.n_trains, params.n_timesteps - 1);
end

function traj = constructTrajectory(params, solution)
    % solution dimensions (train, timestep) values (acceleration -1,0,+1, direction 0-1)
    % initial_pos dimensions (n_trains, 2) values (edge 0-n, position on edge 0-1
    % initial_speed dimensions (n_trains, 1) values (velocity 0-1)
    % trajectory dimensions (train, timestep) values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)

    %% Calculate speed curves
    % Speed is relative to train orientation
    for i_train = 1:params.n_trains
        speeds(i_train, :) = params.initial_speed(i_train) + cumtrapz(solution(i_train, :, 1)) * params.accel;
        # Do not accelerate over maximum
        speeds(i_train, speeds(i_train,:)>params.max_speed) = params.max_speed;
        speeds(i_train,speeds(i_train,:)<-params.max_speed) = -params.max_speed;
    end

    traj(:, 1, :) = params.initial_pos;

    for i_train = 1:params.n_trains
        for timestep = 2:params.n_timesteps-1
            % Start movement in previous state
            traj(i_train, timestep, :) = traj(i_train, timestep-1, :);
            % Add distance to cover this turn
            % Remaining movement is relative to the current edge direction
            remaining_movement = traj(i_train, timestep, 3) * speeds(i_train, timestep);

            abort = 1000;
            while abort > 0
                current_edge_length = params.edge_values(traj(i_train, timestep, 1));
                remaining_backward_length = traj(i_train, timestep, 2) * current_edge_length;
                remaining_forward_length = current_edge_length - remaining_backward_length;

                if (remaining_movement > remaining_forward_length) || (remaining_movement < -remaining_backward_length)
                    % Leave edge forward or backward
                    if remaining_movement > remaining_forward_length
                        traversed_node = params.edge_cols(traj(i_train, timestep, 1));
                        remaining_movement = remaining_movement - remaining_forward_length;
                        node_entrance_direction = traj(i_train, timestep, 3);
                    else
                        traversed_node = params.edge_rows(traj(i_train, timestep, 1));
                        remaining_movement = remaining_movement + remaining_backward_length;
                        node_entrance_direction = -traj(i_train, timestep, 3);
                    end

                    % Find next edge to enter
                    viable_next_edges = find((params.edge_rows == traversed_node) | (params.edge_cols == traversed_node));
                    viable_next_edges = viable_next_edges(viable_next_edges!=traj(i_train, timestep, 1));
                    if length(viable_next_edges) == 0
                        remaining_movement = 0;
                        continue;
                    end

                    % Direction parameter decides next edge
                    next_edge = viable_next_edges(int32(floor(solution(i_train, timestep, 2) * length(viable_next_edges)) + 1));

                    % Move to new edge
                    traj(i_train, timestep, 1) = next_edge;

                    edge_entrance_direction = (params.edge_rows(next_edge) == traversed_node);

                    % Set new edge position and train direction
                    traj(i_train, timestep, 2) = not(edge_entrance_direction);
                    traj(i_train, timestep, 3) = (edge_entrance_direction*2 - 1) * node_entrance_direction;
                    remaining_movement = (edge_entrance_direction*2 - 1) * abs(remaining_movement);

                    abort--;
                else
                    % Adjust position on current edge
                    traj(i_train,timestep, 2) = traj(i_train, timestep, 2) + (remaining_movement / current_edge_length);
                    break;
                end
            end
            assert(abort!=0);
        end
    end
end

function score = objectiveFunction(params, traj)
    penalty = 0;

    %% Separation Penalties
    % For each train pair update the minimum time to collision then skip that time and check again
    for i_train = 1:params.n_trains
        timestep = 1;
        for j_train = i_train+1:params.n_trains
            while timestep < params.n_timesteps
                edge_i = traj(i_train, timestep, 1);
                edge_j = traj(j_train, timestep, 1);
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
