% TODO: needs total refactor
function [position, speeds, start_braking_timestep] = addStop(params, position, speeds, solution, arrival_timestep, departure_time, initial_speed)
    %% Modifies v targets so that train will reach a certain velocity before reaching a position specified by timestep on the old position curve

    start_braking_timestep = findBrakingTimestep(position, speeds, arrival_timestep, params.max_accel);
    [v_target_timesteps, v_target_values] = extractSpeedTargetPoints(params, solution);

    % Consider only points after start of braking
    v_target_relevant_idxs = find(v_target_timesteps >= start_braking_timestep);
    v_target_values = v_target_values(v_target_relevant_idxs);
    v_target_timesteps = v_target_timesteps(v_target_relevant_idxs);

    % Stay stationary during dwell time
    idxs_v_targets_during_stop = (v_target_timesteps >= start_braking_timestep) & (v_target_timesteps <= departure_time);
    v_target_values(idxs_v_targets_during_stop) = 0;

    % Place speed target point of 0 at the start of braking 
    v_target_values(v_target_timesteps == start_braking_timestep) = [];
    v_target_timesteps(v_target_timesteps == start_braking_timestep) = [];
    v_target_timesteps(end+1) = start_braking_timestep;
    v_target_values(end+1) = 0;

    % Recalculate position curve
    if numel(v_target_timesteps) == 1
        v_target = ones(1, params.n_timesteps - start_braking_timestep + 1) * v_target_values(1);
    else
        v_target = interp1(v_target_timesteps, v_target_values, start_braking_timestep:params.n_timesteps, 'previous', 'extrap');
    end
    
    % clf; hold on;
    % plot(position,'DisplayName',"position old");
    % plot(speeds,'DisplayName',"speeds old");

    for i = start_braking_timestep:params.n_timesteps
        if i == 1
            disparity = v_target(i - start_braking_timestep + 1) - initial_speed;
            speeds(i) = initial_speed + sign(disparity) * min(abs(disparity), params.max_accel);
        else
            disparity = v_target(i - start_braking_timestep + 1) - speeds(i-1);
            speeds(i) = speeds(i-1) + sign(disparity) * min(abs(disparity), params.max_accel);
        end
    end

    if start_braking_timestep == 1
        position(start_braking_timestep:params.n_timesteps) = cumsum(speeds(start_braking_timestep:params.n_timesteps));
    else
        position(start_braking_timestep:params.n_timesteps) = position(start_braking_timestep - 1) + cumsum(speeds(start_braking_timestep:params.n_timesteps));
    end

    % plot(position,'DisplayName', "position new");
    % plot(speeds,'DisplayName', "speeds new");
    % scatter(start_braking_timestep, 400,'DisplayName', "start braking timestep");
    % scatter(v_target_timesteps, v_target_values, 'DisplayName',"new speed targets");
    % legend();
end

function start_braking_timestep = findBrakingTimestep(position, speeds, arrival_timestep, max_accel)
    %% Calculate start of braking in order to stop at a position defined by a timestep on the position curve

    approach_direction = sign(speeds(arrival_timestep));
    % Only consider time since last direction change
    first_approach_idx = find(sign(speeds(1:arrival_timestep-1)) == -approach_direction, 1, 'last') + 1;
    if isempty(first_approach_idx)
        first_approach_idx = 1;
    end

    % Select timestep to start braking (will purposefully undershoot due to discretization)
    % Discrete Formula for distance covered under n full braking steps and one remainder braking step
    possible_braking_timesteps = first_approach_idx:arrival_timestep;
    distance_from_stop = abs(position(possible_braking_timesteps) - position(arrival_timestep));
    n = floor(abs(speeds(possible_braking_timesteps)) / max_accel);
    distance_covered_under_braking = n .* (abs(speeds(possible_braking_timesteps)) - 0.5 * (n+1) .* max_accel);
    position_error = distance_covered_under_braking - distance_from_stop;
    
    subset_braking_idxs = find(position_error < 0);

    if isempty(subset_braking_idxs)
        %warning("Failed to undershoot on braking.");
        subset_braking_idxs = 1:numel(position_error);
    end

    [~, best_subset_braking_idx] = min(abs(position_error(subset_braking_idxs)));
    start_braking_timestep = first_approach_idx - 1 + subset_braking_idxs(best_subset_braking_idx);

    if isempty(start_braking_timestep)
        start_braking_timestep = first_approach_idx;
    end
end