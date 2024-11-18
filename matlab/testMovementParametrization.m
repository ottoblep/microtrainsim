params.v_init = 2 * (rand(1,1) - 0.5) * v_max;
params.n_timesteps = 200; % 10s timesteps
params.max_speed = 830; % m/10s = 200km/h
params.max_accel = 46.27; % m/(10s)² = 0-100kmh in 1m
params.n_speed_target_vars = params.n_timesteps / 4; % Support points for acceleration

points = randi([1,params.n_timesteps],1, interp_steps);
vals = rand(1, interp_steps);

params.initial_speed = (rand(1,1) * 2 - 1) * params.max_speed;

constructMovement(params, points, vals, v_init);

function speeds = constructMovement(params, v_target_timesteps, v_target_values, v_init)
    %% Constructs physically possible speed and position curves from target speeds points
    %[v_target_timesteps, v_target_sorted_idxs] = sort(v_target_timesteps);
    %v_target_values = v_target_values(v_target_sorted_idxs);

    % Fill timesteps
    v_target_timesteps(1) = 1;
    [v_target_timesteps , unique_idxs, ~] = unique(v_target_timesteps);
    v_target_values = (2 * v_target_values(unique_idxs) - 1) * params.max_speed;

    % Interpolate
    v_target = interp1(v_target_timesteps, v_target_values, 1:params.n_timesteps, 'previous', 'extrap');

    speeds = zeros(1, params.n_timesteps);
    speeds(1) = params.initial_speed;
    for i = 2:params.n_timesteps
        diff = v_target(i) - speeds(i-1);
        speeds(i) = speeds(i-1) + sign(diff) * min(abs(diff), params.max_accel);
    end

    clf; hold on;
    scatter(v_target_timesteps, v_target_values);
    plot(v_target);
    plot(speeds);
end