function [v_targets, start_braking_timestep, end_braking_timestep] = addBraking(params, global_speeds, v_targets, edge_transition, immediate_departure, hold_until_timestep, braking_goal_speed)
    %% Returns a subset of v_targets modified so that train will reach a certain velocity before reaching the new edge specified in edge_transition
    % Moves speed to braking_goal speed regardless of direction changes 

    assert(abs(edge_transition.speed) >= abs(braking_goal_speed));
    % global_speeds ranges from start of simulation until edge_transition.timestep
    assert(numel(global_speeds) == edge_transition.timestep);

    % Ensure directionality
    braking_goal_speed = sign(edge_transition.node_traversal_direction) * abs(braking_goal_speed);

    v_target_timesteps = v_targets(:, 1);
    v_target_values = v_targets(:, 2);

    [start_braking_timestep end_braking_timestep] = findBrakingTimestep(params, global_speeds, edge_transition, braking_goal_speed);

    % Remove other v targets during braking 
    if not(immediate_departure) & (hold_until_timestep > end_braking_timestep)
        end_braking_timestep = hold_until_timestep;
    end

    idxs_v_targets_during_braking = find((v_target_timesteps >= start_braking_timestep) & (v_target_timesteps <= end_braking_timestep));
    v_target_values(idxs_v_targets_during_braking) = [];
    v_target_timesteps(idxs_v_targets_during_braking) = [];

    % Place speed target point of 0 at the start of braking or not accelerating
    v_target_timesteps(end+1) = start_braking_timestep;
    v_target_values(end+1) = braking_goal_speed;

    v_targets = [v_target_timesteps v_target_values];
end

function [start_braking_timestep, end_braking_timestep] = findBrakingTimestep(params, global_speeds, edge_transition, braking_goal_speed)
    %% Find start_braking_timestep that undershoots the position at edge_transition.timestep the littlest while braking with max_accel starting at start_braking_timestep
    
    global_trajectory = zeros(1, edge_transition.timestep);
    global_trajectory(2:edge_transition.timestep) = cumsum(global_speeds(1:edge_transition.timestep - 1));
    target_position = global_trajectory(edge_transition.timestep) - edge_transition.node_traversal_direction * edge_transition.extra_movement;

    % Select timestep to start braking 
    % Will purposefully undershoot on position due to discretization
    % Candidate timesteps are the steps where braking actually takes effect 
    % start_braking_timestep (speed target point) is 1 timestep earlier due to propagation delay
    candidate_timesteps = 2:edge_transition.timestep;
    abs_diff_to_target_speed = abs(global_speeds(candidate_timesteps) - braking_goal_speed);
    required_braking_time = ceil(abs_diff_to_target_speed / params.max_accel);
    goal_speed_reached_at = candidate_timesteps + required_braking_time';
    accel_direction = sign(braking_goal_speed - global_speeds(candidate_timesteps));
    % We need a general integral here including zero crossings that can happen
    dist_covered_while_braking = global_speeds(candidate_timesteps) .* required_braking_time ...
                                 + accel_direction .* (0.5 * (required_braking_time + 1) .* required_braking_time * params.max_accel);
    % Overshot = negative dist remaining
    dist_remaining_after_braking = target_position - global_trajectory(candidate_timesteps) - dist_covered_while_braking';
    
    % Find latest possible braking point
    % Speed goal must be reached the step before the transition (p(n) = n(n-1) + v(n-1))
    best_candidate_braking_idx = find(dist_remaining_after_braking >= 0 & goal_speed_reached_at <= edge_transition.timestep, 1, 'last');

    start_braking_timestep = candidate_timesteps(best_candidate_braking_idx) - 1;

    if (isempty(start_braking_timestep) || start_braking_timestep > edge_transition.timestep || start_braking_timestep < 1)
        error("Failed to find braking timestep.");
    end

    end_braking_timestep = max(goal_speed_reached_at(best_candidate_braking_idx), edge_transition.timestep);
end
