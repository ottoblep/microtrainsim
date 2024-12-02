[network, params] = generateEnvironment("crossover");

single_train_solution = rand(1, 2 * params.n_v_target_vars + params.n_switch_vars);

[traj, events, n_fullfilled_stops] = constructTrajectory(network, params, params.initial_positions(1, :), params.initial_speeds(1), params.planned_stops(params.planned_stops(:,1)==1, 2:3));