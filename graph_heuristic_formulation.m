params.infra = [ 0 50 10 0;
                 0 0 60 0;
                 0 0 0 40;
                 0 0 0 0];
[params.edge_rows, params.edge_cols, params.edge_values] = find(params.infra);
params.n_trains = 3;
params.initial_pos = [2, 0.5; 3, 0.5; 1,0.1;];
params.initial_speed = [0.2, 0.3, 0.4];
params.n_timesteps = 100;
params.min_separation = 5;
params.max_speed = 1;
params.accel = 0.1;

function sol = randomSolution(params)
    % Solution Array dimensions (train, timestep)
    % Solution Array values (accelerate -1,0,1, direction 0-1)
    sol = rand(params.n_trains, params.n_timesteps - 1, 2);
    sol(:,:,1) = randi([-1,1], params.n_trains, params.n_timesteps - 1);
end

function traj = constructTrajectory(params, solution)
    % solution dimensions (train, timestep) values (target_speed 0-1, direction 0-1)
    % initial_pos dimensions (n_trains, 2) values (edge 0-n, position_on_edge 0-1
    % initial_speed dimensions (n_trains, 1) values (velocity 0-1)
    % trajectory dimensions (train, timestep) values (edge 0-n, position_on_edge 0-1)
    traj(:, 1, :) = params.initial_pos;
    speeds(:, 1) = params.initial_speed;

    for timestep = 2:params.n_timesteps-1
        for i_train = 1:params.n_trains
            %% Calculate speeds
            new_speed = speeds(i_train, timestep-1) + params.accel * solution(i_train, timestep, 1);
            if abs(new_speed) < params.max_speed
                speeds(i_train, timestep) = new_speed;
            else
                speeds(i_train, timestep) = speeds(i_train, timestep-1);
            end

            %% Move trains
            % Start movement at previous edge
            traj(i_train, timestep, 1) = traj(i_train, timestep-1, 1);
            remaining_movement = speeds(i_train, timestep);

            abort = 10000;
            while abort > 0
                if abort < 10
                   disp("oh oh");
                end
                current_edge_length = params.edge_values(traj(i_train, timestep, 1));

                remaining_forward_length = (1 - traj(i_train, timestep, 2)) * current_edge_length;
                remaining_backward_length = traj(i_train, timestep, 2) * current_edge_length;

                if abs(remaining_movement) < 1e-10
                    remaining_movement = 0;
                end

                if (remaining_movement > remaining_forward_length) || (remaining_movement < -remaining_backward_length)
                    % Leave edge forward or backward
                    if remaining_movement > remaining_forward_length
                        traversed_node = params.edge_rows(traj(i_train, timestep, 1));
                        remaining_movement = remaining_movement - remaining_forward_length;
                    else
                        traversed_node = params.edge_cols(traj(i_train, timestep, 1));
                        remaining_movement = remaining_movement + remaining_backward_length;
                    end

                    viable_next_edges = find((params.edge_rows == traversed_node) | (params.edge_cols == traversed_node));
                    viable_next_edges = viable_next_edges(viable_next_edges!=traj(i_train, timestep, 1));
                    if length(viable_next_edges) == 0
                        traj(i_train, timestep, 2) = (sign(remaining_movement) + 1) / 2;
                        remaining_movement = 0;
                        speeds(i_train, timestep) = 0;
                        continue;
                    end

                    % Direction parameter decides next edge
                    next_edge = viable_next_edges(int32(floor(solution(i_train, timestep, 2) * length(viable_next_edges)) + 1));

                    % Move to new edge
                    traj(i_train, timestep, 1) = next_edge;

                    % Set new edge position and movement direction 
                    if params.edge_rows(next_edge) == traversed_node
                        traj(i_train, timestep, 2) = 0;
                        remaining_movement = abs(remaining_movement);
                        speeds(i_train, timestep) = abs(speeds(i_train, timestep));
                    else
                        traj(i_train, timestep, 2) = 1;
                        remaining_movement = -abs(remaining_movement);
                        speeds(i_train, timestep) = -abs(speeds(i_train, timestep));
                    end

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

function score = objectiveFunction(solution, params)
    penalty = 0;

    %% Separation Penalties
    for timestep = 1:params.n_timesteps
        for i_train = 1:params.n_trains
            distance = boundedDistance(params, i_train, solution, timestep);

            % Penalize smaller than safe distances
            if distance >= params.min_separation
                continue;
            else
                penalty = penalty + params.min_separation - distance;
            end
        end
    end

    %% Objective evaluation
    score = -penalty;
end

function dist = boundedDistance(params, edge1, position1, edge2, position2)
    hgraph = params.infra;
    old_size = size(params.infra)(1);

    if edge1 == edge2
        dist = abs((position2-position1) * params.infra(params.edge_rows(edge1),params.edge_cols(edge1)));
    % else
        % TODO: Insert helper nodes
        % TODO: Bounded BFS
    end
end

tic
sol = randomSolution(params);
traj = constructTrajectory(params, sol);
toc
