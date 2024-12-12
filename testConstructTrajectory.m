[network, params] = generateEnvironment("crossover");

for i = 1:1000
    single_train_solution = rand(1, 2 * params.n_v_target_vars + params.n_switch_vars);
    [traj, events] = constructTrajectory(network, params, single_train_solution, params.initial_positions(1, :), params.initial_speeds(1), params.planned_stops((params.planned_stops(:,1)==1), 2:3));

    % Parameter bounds
    assert(not(any(any(isnan(traj)))));
    assert(all(traj(:,2) >= 0) && all(traj(:,2) <= 1));
    % Speed limits
    assert(all(abs(traj(:,4))' <= network.speed_limits(traj(:,1)) + 1e-14));
end
