params.n_timesteps = 200; % 10s timesteps
params.max_speed = 830; % m/10s = 200km/h
params.max_accel = 46.27; % m/(10s)² = 0-100kmh in 1m
params.n_speed_target_vars = params.n_timesteps / 4; % Support points for acceleration
params.n_switch_vars = 20; % unused dummy
initial_position = 3;
initial_speed = (rand(1,1) * 2 - 1) * params.max_speed;
solution = [randi([1,params.n_timesteps],1, params.n_speed_target_vars), rand(1, params.n_speed_target_vars), rand(1, params.n_switch_vars)]; % dummy
arrival_timestep = 50;
departure_timestep = 100;

speeds = constructMovement(params, solution, initial_speed);
position = cumsum(speeds);
solution_new = addStop(params, position, speeds, solution, 50, 100, initial_speed, initial_position, 1);
speeds_new = constructMovement(params, solution_new, initial_speed);
position_new = cumsum(speeds_new);

clf; hold on;
scatter(arrival_timestep, position(arrival_timestep));
plot(speeds);
plot(speeds_new);
plot(position);
plot(position_new);

function solution = addStop(params, position, speeds, solution, arrival_timestep, departure_time, initial_speed, initial_position, overshoot)
    %% Modifies solution to stop around a certain position defined by a timestep on the old position curve

    approach_direction = sign(speeds(arrival_timestep));
    % Only consider time since last direction change
    first_approach_idx = find(sign(speeds(1:arrival_timestep)) ~= approach_direction, 1, 'last') + 1;
    if isempty(first_approach_idx)
        first_approach_idx = 1;
    end

    % Select timestep to start braking (will over or undershoot due to discretization)
    % Discrete Formula for distance covered under n full braking steps and
    % one remainder braking step
    possible_braking_timesteps = first_approach_idx:arrival_timestep;
    distance_from_stop = abs(position(possible_braking_timesteps) - position(arrival_timestep));
    n = floor(abs(speeds(possible_braking_timesteps)) / params.max_accel);
    distance_covered_under_braking = n .* (abs(speeds(possible_braking_timesteps)) - 0.5 * (n+1) .* params.max_accel);
    [~, best_braking_idx] = min(abs(distance_covered_under_braking - distance_from_stop));
    start_braking_timestep = first_approach_idx + best_braking_idx;
    if isempty(start_braking_timestep)
        start_braking_timestep = first_approach_idx;
    end
    
    % Adjust acceleration points
    v_target_timesteps = solution(1:params.n_speed_target_vars);
    v_target_values = solution(params.n_speed_target_vars + 1:2 * params.n_speed_target_vars);
    [v_target_timesteps, v_target_sorted_idxs] = sort(v_target_timesteps);
    v_target_values = v_target_values(v_target_sorted_idxs);

    % Stay stationary 
    idxs_v_targets_during_stop = (v_target_timesteps >= start_braking_timestep) & (v_target_timesteps <= departure_time);
    v_target_values(idxs_v_targets_during_stop) = 0.5;
    % Start braking instantly
    idx_v_target_to_shift = find(v_target_timesteps > start_braking_timestep, 1, 'first');
    v_target_timesteps(idx_v_target_to_shift) = start_braking_timestep;
    v_target_values(idx_v_target_to_shift) = 0.5;

    solution(1:params.n_speed_target_vars) = v_target_timesteps;
    solution(params.n_speed_target_vars + 1:2 * params.n_speed_target_vars) = v_target_values;
end

function speeds = constructMovement(params, solution, initial_speed)
    %% Constructs physically possible speed and position curves from target speeds points
    v_target_timesteps = solution(1:params.n_speed_target_vars);
    v_target_values = solution(params.n_speed_target_vars + 1:2 * params.n_speed_target_vars);

    [~, first_v_target_timestep_idx] = min(v_target_timesteps);
    v_target_timesteps(first_v_target_timestep_idx) = 1;
    [v_target_timesteps , unique_idxs, ~] = unique(v_target_timesteps);
    v_target_values = (2 * v_target_values(unique_idxs) - 1) * params.max_speed;
    v_target = interp1(v_target_timesteps, v_target_values, 1:params.n_timesteps, 'previous');

    speeds = zeros(1, params.n_timesteps);
    speeds(1) = initial_speed;
    for i = 2:params.n_timesteps
        diff = v_target(i) - speeds(i-1);
        speeds(i) = speeds(i-1) + sign(diff) * min(abs(diff), params.max_accel);
    end

    % clf; hold on;
    % scatter(v_target_timesteps, v_target_values);
    % plot(v_target);
    % plot(speeds);
end