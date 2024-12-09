function [v_targets, start_braking_timestep, end_braking_timestep] = addBraking(params, global_speeds, v_targets, edge_transition, immediate_departure, hold_until_timestep, braking_goal_speed)
    %% Returns a subset of v_targets modified so that train will reach a certain velocity before reaching the new edge specified in edge_transition
    assert(abs(edge_transition.speed) > abs(braking_goal_speed));
    % global_speeds ranges from start of simulation until edge_transition.timestep
    assert(numel(global_speeds) == edge_transition.timestep);

    v_target_timesteps = v_targets(:, 1);
    v_target_values = v_targets(:, 2);

    [start_braking_timestep end_braking_timestep] = findBrakingTimestep(params, global_speeds, edge_transition, braking_goal_speed);

    % Remove other v targets during braking 
    if not(immediate_departure) && (hold_until_timestep > end_braking_timestep)
        end_braking_timestep = hold_until_timestep;
    end

    idxs_v_targets_during_braking = find((v_target_timesteps >= start_braking_timestep) & (v_target_timesteps < end_braking_timestep));
    v_target_values(idxs_v_targets_during_braking) = [];
    v_target_timesteps(idxs_v_targets_during_braking) = [];

    % Place speed target point of 0 at the start of braking 
    v_target_timesteps(end+1) = start_braking_timestep;
    v_target_values(end+1) = edge_transition.node_traversal_direction * braking_goal_speed;

    v_targets = [v_target_timesteps v_target_values];
end

function [start_braking_timestep, end_braking_timestep] = findBrakingTimestep(params, global_speeds, edge_transition, braking_goal_speed)
    %% Find start_braking_timestep that undershoots the position at edge_transition.timestep the littlest while braking with max_accel starting at start_braking_timestep

    global_trajectory = zeros(1, edge_transition.timestep);
    global_trajectory(2:edge_transition.timestep) = cumsum(global_speeds(1:edge_transition.timestep - 1));

    % Only consider time since last direction change
    % first_approach_idx = find(sign(global_speeds) == -edge_transition.node_traversal_direction, 1, 'last') - 1;
    % if isempty(first_approach_idx)
    %     first_approach_idx = 1;
    % end
    first_approach_idx = 1;

    % Select timestep to start braking (will purposefully undershoot on position due to discretization)
    candidate_timesteps = first_approach_idx:edge_transition.timestep - 1;
    target_position = global_trajectory(edge_transition.timestep) - edge_transition.node_traversal_direction * edge_transition.extra_movement;
    abs_v_diff = abs(global_speeds(candidate_timesteps)) - abs(braking_goal_speed);
    required_braking_time = ceil(abs_v_diff / params.max_accel);
    speed_reached_at = candidate_timesteps + required_braking_time';
    dist_covered_while_braking = required_braking_time * braking_goal_speed + (0.5 * (required_braking_time - 1) * params.max_accel) + rem(abs_v_diff, params.max_accel);
    dist_remaining_after_braking = abs(global_trajectory(candidate_timesteps) - target_position) - dist_covered_while_braking';
    
    % Undershoot as little as possible
    % Speed goal must be reached the step before the transition (p(n) = n(n-1) + v(n-1))
    subset_braking_idxs = find(dist_remaining_after_braking > 0 & speed_reached_at < edge_transition.timestep - 1);

    if isempty(subset_braking_idxs)
        error("Failed to undershoot on braking.");
        % subset_braking_idxs = 1:numel(dist_remaining_after_braking);
    end

    [~, best_subset_braking_idx] = min(abs(dist_remaining_after_braking(subset_braking_idxs)));
    start_braking_timestep = candidate_timesteps(subset_braking_idxs(best_subset_braking_idx));
    end_braking_timestep = start_braking_timestep + required_braking_time(subset_braking_idxs(best_subset_braking_idx));

    if (isempty(start_braking_timestep) || start_braking_timestep > edge_transition.timestep || start_braking_timestep < 1)
        error("Failed to find braking timestep.");
        % start_braking_timestep = first_approach_idx;
    end
end
