function [v_targets start_braking_timestep] = addBraking(params, global_speeds, v_targets, previous_arrival_timestep, immediate_departure, hold_until_timestep, braking_goal_speed)
    %% Modifies v targets so that train will reach a certain velocity before reaching a position specified by timestep on the old position curve
    v_target_timesteps = v_targets(:, 1);
    v_target_values = v_targets(:, 2);

    start_braking_timestep = findBrakingTimestep(params, global_speeds, previous_arrival_timestep, braking_goal_speed);

    % Consider only points after start of braking
    v_target_relevant_idxs = find(v_target_timesteps >= start_braking_timestep);
    v_target_values = v_target_values(v_target_relevant_idxs);
    v_target_timesteps = v_target_timesteps(v_target_relevant_idxs);

    % Stay stationary until departure
    % Dwell time only possible when stationary
    if not(immediate_departure)
        idxs_v_targets_during_stop = (v_target_timesteps >= start_braking_timestep) & (v_target_timesteps < hold_until_timestep);
        v_target_values(idxs_v_targets_during_stop) = braking_goal_speed;
    end

    % Place speed target point of 0 at the start of braking 
    v_target_values(v_target_timesteps == start_braking_timestep) = [];
    v_target_timesteps(v_target_timesteps == start_braking_timestep) = [];
    v_target_timesteps(end+1) = start_braking_timestep;
    v_target_values(end+1) = braking_goal_speed;

    v_targets = [v_target_timesteps v_target_values];
end

function start_braking_timestep = findBrakingTimestep(params, global_speeds, previous_arrival_timestep, braking_goal_speed)
    %% Find start_braking_timestep that undershoots the position at previous_arrival_timestep the littlest while braking with max_accel starting at start_braking_timestep

    % global_speeds speeds from start of simulation until previous_arrival_timestep
    assert(numel(global_speeds) == previous_arrival_timestep);

    global_trajectory = zeros(1,previous_arrival_timestep);
    global_trajectory(2:previous_arrival_timestep) = cumsum(global_speeds(1:previous_arrival_timestep - 1));

    approach_direction = sign(global_speeds(end));
    % Only consider time since last direction change
    first_approach_idx = find(sign(global_speeds) == -approach_direction, 1, 'last') + 1;
    if isempty(first_approach_idx)
        first_approach_idx = 1;
    end

    % Select timestep to start braking (will purposefully undershoot on position due to discretization)
    candidate_timesteps = first_approach_idx:previous_arrival_timestep - 1;
    target_position = global_trajectory(previous_arrival_timestep);
    abs_v_diff = abs(global_speeds(candidate_timesteps) - braking_goal_speed);
    required_braking_time = ceil(abs_v_diff / params.max_accel);
    distance_covered_while_braking = required_braking_time * braking_goal_speed + (0.5 * (required_braking_time - 1) * params.max_accel) + rem(abs_v_diff, params.max_accel);
    position_error = abs(global_trajectory(candidate_timesteps) - target_position) - distance_covered_while_braking';
    
    subset_braking_idxs = find(position_error <= 0);

    if isempty(subset_braking_idxs)
        warning("Failed to undershoot on braking.");
        subset_braking_idxs = 1:numel(position_error);
    end

    [~, best_subset_braking_idx] = min(abs(position_error(subset_braking_idxs)));
    start_braking_timestep = first_approach_idx - 1 + candidate_timesteps(subset_braking_idxs(best_subset_braking_idx));

    if isempty(start_braking_timestep)
        warning("Failed to find braking timestep.");
        start_braking_timestep = first_approach_idx;
    end
end
