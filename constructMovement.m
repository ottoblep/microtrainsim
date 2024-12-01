function speeds = constructMovement(v_targets, start_timestep, end_timestep, initial_speed, max_accel)
    %% Constructs physically possible speed curve from target speeds points

    v_target = interp1(v_target_timesteps, v_target_values, start_timestep:end_timestep, 'previous', 'extrap');
    v_target(1) = v_target(1) - (start_timestep - 1);

    speeds = zeros(1, end_timestep - start_timestep + 1);
    speeds(1) = inital_speed;
    for i = 2:(end_timestep - start_timestep + 1)
        disparity = v_target(i-1) - speeds(i-1);
        speeds(i) = speeds(i-1) + sign(disparity) * min(abs(disparity), max_accel);
    end

    % clf; hold on;
    % scatter(v_target_timesteps, v_target_values);
    % plot(v_target);
    % plot(speeds);
end