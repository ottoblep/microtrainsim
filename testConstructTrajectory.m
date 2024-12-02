[network, params] = generateEnvironment("crossover");

single_train_solution = rand(1, 2 * params.n_v_target_vars + params.n_switch_vars);

[traj, events] = constructTrajectory(network, params, single_train_solution, params.initial_positions(1, :), params.initial_speeds(1), params.planned_stops((params.planned_stops(:,1)==1), 2:3));