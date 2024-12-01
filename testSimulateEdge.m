v_targets = [[0.1 0.3]; [0.4 -0.5]; [0.6 0.8]; [0.9, 0.2]];
[network, params] = generateEnvironment("crossover");
v_targets(:,1) = v_targets(:,1) * params.n_timesteps;
v_targets(:,2) = v_targets(:,2) * params.max_speed;
initial_edge_state = [35 3 0.2 -1 20];

[edge_transition edge_trajectory] = simulateEdge(network, params, initial_edge_state, v_targets);