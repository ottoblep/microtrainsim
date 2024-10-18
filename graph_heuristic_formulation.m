% Generate Network
adj = random_planar_graph(10);
network.adjacency_matrix = triu((adj + adj') * 1000, 1);
%network.adjacency_matrix = [ 0 100;
%                             0 0 ];
connection_indexes = find(network.adjacency_matrix);
network.adjacency_matrix(connection_indexes) = network.adjacency_matrix(connection_indexes) .* (randi([1,3],size(connection_indexes)));
clear adj connection_indexes;
[network.edge_rows, network.edge_cols, network.edge_values] = find(network.adjacency_matrix);
network.adjacent_edge_list = {};
for node = 1:size(network.adjacency_matrix,1)
    network.adjacent_edge_list{node} = find((network.edge_rows == node) | (network.edge_cols == node));
end

%% All shortest path pairs
tmp_adj = network.adjacency_matrix;
tmp_adj(tmp_adj==0) = Inf;
network.all_shortest_paths = distances(graph(tmp_adj,'upper'), 'Method', 'positive');
clear tmp_adj;

%% Parameters
params.n_timesteps = 8640; % 10s timesteps for one whole day
params.n_trains = 3;
params.min_separation = 100; % m
params.max_speed = 1.11; % m/10s 200km/h
params.max_accel = 46.27; % m/(10s)² 0-100kmh in 1m
params.max_changeover_time = 1440; % 4hrs
params.train_capacity = 400; % 400 Passengers

params.demand_matrix = randi(1000, size(network.adjacency_matrix,1));
params.demand_matrix(logical(eye(size(params.demand_matrix)))) = 0;
%params.demand_matrix = [ 0 20;
%                         0 0];
params.initial_positions = randi([1,length(network.edge_values)], params.n_trains, 3);
params.initial_positions(:,2) = rand(params.n_trains,1);
params.initial_positions(:,3) = randi([0,1],params.n_trains,1) * 2 - 1;
params.initial_speeds = (rand(params.n_trains,1) * 2 - 1) * params.max_speed;

%% Solution Construction

function [traj_set, event_set] = constructTrajectorySet(network, solution, initial_positions, initial_speeds, max_accel, max_speed)
    %% Constructs a set of train trajectories on the graph
    % solution dimensions (n_trains, timestep * 2)
    % solution values (acceleration 0-1, direction 0-1)
    % initial position dimensions (n_trains, 3)
    % initial position values (edge 0-n, position on edge 0-1)
    % initial speed dimensions (n_trains)
    % initial speed values (speed: -max_speed-max_speed)
    % traj_set dimensions (n_trains, 3, timestep)
    % traj_set values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)
    % event_set dimensions (3, n))
    % event_set values (presence_at_node, timestep, train)
    n_trains = size(solution,1);
    event_set = [];
    for i_train = 1:n_trains
        [traj_set(i_train, :, :), new_events] = constructTrajectory(network, solution(i_train,:), initial_positions(i_train, :), initial_speeds(i_train), max_accel, max_speed);
        new_events(3,:) = i_train;
        event_set = cat(2, event_set, new_events);
    end
end

function [traj, events] = constructTrajectory(network, solution, initial_position, initial_speed, max_accel, max_speed)
    %% Constructs a single train trajectory on the graph
    % solution dimensions (timestep * 2)
    % solution values (acceleration 0-1, direction 0-1)
    % initial position dimensions (3)
    % initial position values (edge 0-n, position on edge 0-1)
    % initial speed dimensions (1)
    % initial speed values (speed: -max_speed-max_speed)
    % trajectory dimensions (3, timestep)
    % trajectory values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)
    % events dimensions (2, n)
    % events values (presence_at_node, timestep)
    assert(~any(solution>1));
    assert(~any(solution<0));

    n_timesteps = size(solution,2) / 2;

    % Calculate speed and position curves relative to train orientation
    acceleration = solution(1:n_timesteps);
    speeds = initial_speed + cumtrapz(((acceleration * 2) - 1) * max_accel);
    % Do not accelerate over maximum
    speeds(speeds>max_speed) = max_speed;
    speeds(speeds<-max_speed) = -max_speed;
    position = cumtrapz(speeds);

    events = [];
    traj(1, :) = initial_position(1);
    traj(2, :) = initial_position(2);
    traj(3, :) = initial_position(3);

    pivot_timestep = 1;
    while pivot_timestep < n_timesteps
        current_edge_length = network.edge_values(traj(1, pivot_timestep));
        remaining_backward_length = traj(2, pivot_timestep) * current_edge_length;
        remaining_forward_length = current_edge_length - remaining_backward_length;
        base_position = position(pivot_timestep);

        % Find next edge change
        edge_trajectory = traj(3, pivot_timestep) * (position(pivot_timestep:end) - base_position);
        next_pivot_timestep = (pivot_timestep - 1) + find((edge_trajectory > remaining_forward_length | edge_trajectory < -remaining_backward_length), 1);
        edge_exit_point = any((edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) > remaining_forward_length));
        if isempty(next_pivot_timestep)
            next_pivot_timestep = n_timesteps;
        end
        assert(next_pivot_timestep <= n_timesteps);

        % Set trajectory on the current edge
        traj(1, pivot_timestep:next_pivot_timestep-1) = traj(1, pivot_timestep);
        traj(2, pivot_timestep:next_pivot_timestep-1) = traj(2, pivot_timestep) + edge_trajectory(1:next_pivot_timestep - (pivot_timestep-1) - 1) / current_edge_length;
        traj(3, pivot_timestep:next_pivot_timestep-1) = traj(3, pivot_timestep);

        % Leave current edge
        % node_traversal_direction = edge_exit_point XNOR old_train_orientation
        if edge_exit_point
            traversed_node = network.edge_cols(traj(1, pivot_timestep));
            node_traversal_direction = traj(3, pivot_timestep);
            extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) - remaining_forward_length;
        else
            traversed_node = network.edge_rows(traj(1, pivot_timestep));
            node_traversal_direction = -traj(3, pivot_timestep);
            extra_movement = edge_trajectory(next_pivot_timestep - (pivot_timestep - 1)) + remaining_backward_length;
        end

        events(:, end+1) = [traversed_node, next_pivot_timestep];

        % Decide on next edge
        viable_next_edges = network.adjacent_edge_list{traversed_node};
        viable_next_edges = viable_next_edges(viable_next_edges~=traj(1, pivot_timestep));

        if isempty(viable_next_edges)
            % Stay stationary until trajectory comes back around
            % Use the old base position and forw/backw lengths
            deadend_edge_trajectory = traj(3, pivot_timestep) * (position(next_pivot_timestep:end) - base_position);
            if edge_exit_point
                deadend_next_pivot_timestep = (next_pivot_timestep - 1) + find(deadend_edge_trajectory < remaining_forward_length, 1);
            else
                deadend_next_pivot_timestep = (next_pivot_timestep - 1) + find(deadend_edge_trajectory > -remaining_backward_length, 1);
            end
            if isempty(deadend_next_pivot_timestep)
                deadend_next_pivot_timestep = n_timesteps;
            end
            assert(deadend_next_pivot_timestep <= n_timesteps);

            % Set position and train direction as stationary
            traj(1, next_pivot_timestep:deadend_next_pivot_timestep) = traj(1, pivot_timestep);
            traj(2, next_pivot_timestep:deadend_next_pivot_timestep) = double(edge_exit_point);
            traj(3, next_pivot_timestep:deadend_next_pivot_timestep) = traj(3, pivot_timestep);

            pivot_timestep = deadend_next_pivot_timestep;
        else
            next_edge_selection = 1 + round(solution(n_timesteps + pivot_timestep) * (length(viable_next_edges) - 1));
            next_edge = viable_next_edges(next_edge_selection);

            edge_entrance_point = (network.edge_cols(next_edge) == traversed_node);
            % New Edge
            traj(1, next_pivot_timestep) = next_edge;
            % New Position on Edge
            traj(2, next_pivot_timestep) = edge_entrance_point + (-1) * (edge_entrance_point*2 - 1) * abs(extra_movement / network.edge_values(next_edge));
            % New Orientation on Edge
            % new_train_orientation      =         edge_entrance_point      XOR node_traversal_direction    ( XNOR is multiplication )
            traj(3, next_pivot_timestep) = (-1) * (edge_entrance_point*2 - 1) * node_traversal_direction;

            pivot_timestep = next_pivot_timestep;
        end
    end

    traj(:, n_timesteps) = traj(:, n_timesteps-1);
end

%% Objective Evalutation

function collision_score = collisionPenalties(network, traj_set, min_separation, max_speed)
    %% Evaluates a set of train trajectories
    % trajectory_set dimensions (n_trains, 3, timestep)
    % trajectory_set values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)
    penalty = 0;

    % Separation Penalties
    % For each train pair update the minimum time to collision then skip that time and check again
    n_trains = size(traj_set, 1);
    n_timesteps = size(traj_set, 3);
    parfor i_train = 1:n_trains
        for j_train = i_train+1:n_trains
            timestep = 1;
            while timestep < n_timesteps
                distance = trainDistance(network, traj_set, i_train, j_train, timestep);

                if distance > min_separation
                    guaranteed_safe_time = int32(floor((distance - min_separation) / (2 * max_speed))) + 1;
                else
                    % Exponential penalty for closeness beyond the minimum separation
                    penalty = penalty + min(10000, (min_separation/distance - 1));
                    guaranteed_safe_time = 1;
                end

                timestep = timestep + guaranteed_safe_time;
            end
        end
    end

    collision_score = -penalty;
end

function [demand_score, transfer_graph_digraph, flow_solution] = demandSatisfaction(network, event_set, demand_matrix, max_changeover_time, train_capacity)
    %% Evaluates satisfied demand for a given set of arrival/departure timesteps 
    % event_set dimenstions (3, n))
    % event_set values (presence_at_node, timestep, train)
    assert(all(size(demand_matrix) == size(network.adjacency_matrix)));

    transfer_graph = constructTransferGraph(network, event_set, max_changeover_time, train_capacity);
    transfer_graph_digraph = digraph(transfer_graph);
    transfer_graph_edge_idxs = find(transfer_graph);

    n_nodes = size(transfer_graph,1); % n_stations + n_stops + n_stations
    n_edges = numel(transfer_graph_edge_idxs);
    n_stations = size(network.adjacency_matrix,1);
    % Only nondiagonal entries in the demand matrix
    n_demands = n_stations^2 - n_stations;
    n_decision_vars = (n_edges + 1)* n_demands;

    %% Splittable Multi-commodity maximum flow problem https://en.wikipedia.org/wiki/Multi-commodity_flow_problem
    %  NP-Complete if discrete but polynomial if real
    %  no loops, directed, positive weights

    % Decision Variables:
    % - edge flows, real positive, per demand, per edge
    % - carried demand, real positive, per demand relation 

    % Objective:
    % - maximize sum of carried demand
    f = cat(1, zeros(n_edges * n_demands, 1), -1 * ones(n_demands, 1));

    % Constraints:
    % - vertex flow conservation constraint, equality, per demand relation, per node

    %  matrix size is demands approx stations^4 * journeys^2
    %  sparse matrix contains 2 * (stations^2 - stations) * (journeys + 1) nonzero elements
    %                           | Edge_flows D1 | Edge_flows D2 | carried demand | 
    % Demand 1          Sources |               |               | -1  0  | 
    %           Nodes   Stops   |      A        |       0       | 0 0  | 
    %                   Sinks   |               |               | 1  0  | =   0
    % Demand 2          Sources |               |               | 0  -1  | 
    %           Nodes   Stops   |      0        |       B       | 0  0 | 
    %                   Sinks   |               |               | 0  1  | 
    Aeq = sparse(n_demands * n_nodes,n_decision_vars);
    beq = sparse(n_demands * n_nodes, 1);

    % General flow constraint template in the network (A)
    Aeq_single_flow = sparse(n_nodes, n_edges);
    for i_node = 1:n_nodes
        out_edges_idxs = outedges(transfer_graph_digraph, i_node);
        in_edges_idxs = inedges(transfer_graph_digraph, i_node);
        Aeq_single_flow(i_node, in_edges_idxs) = 1; % inflows 
        Aeq_single_flow(i_node, out_edges_idxs) = -1; % outflows 
    end
    % source and sink nodes
    assert(all(all(Aeq_single_flow(1:n_stations, :) <= 0)));
    assert(all(all(Aeq_single_flow(end - n_stations + 1:end, :) >= 0)));

    % Copy flow constraints per demand
    k_demand = 1;
    for j_demand_source = 1:n_stations
        for i_demand_destination = 1:n_stations
            if i_demand_destination ==  j_demand_source
                continue
            end

            % copy flow constraint matrix (A,B)
            Aeq((k_demand-1) * n_nodes + 1:k_demand * n_nodes, (k_demand-1) * n_edges + 1:k_demand * n_edges) = Aeq_single_flow;
            % carried demand source node
            Aeq((k_demand-1) * n_nodes + j_demand_source, n_demands * n_edges + k_demand) = 1;
            % carried demand sink node
            Aeq(k_demand * n_nodes - n_stations + i_demand_destination, n_demands * n_edges + k_demand) = -1;

            k_demand = k_demand + 1;
        end
    end
    % column sum is always zero 
    assert(all(sum(Aeq) == zeros(1, size(Aeq,2))));

    % - edge capacity constraint, inequality, per edge
    % - available demand constraint, inequality, per demand

    %                 | Edges * Demands | Demands |
    % Edges * Demands |   I             |   0     | <= A
    % Demands         |   0             |   I     | <= B
    A = speye(n_decision_vars);
    b = zeros(n_decision_vars, 1);
    b(1:n_edges*n_demands) = repmat(transfer_graph(transfer_graph_edge_idxs), n_demands, 1);
    demands_transposed = transpose(demand_matrix);
    b(end - n_demands + 1:end) = demands_transposed(~eye(size(demand_matrix)));
    lb = sparse(n_decision_vars, 1);
    options = optimoptions('linprog','Display','none');
    % Run Solver
    [flow_solution, obj_val, ~, output] = linprog(f, A, b, Aeq, beq, lb, [], options);
    demand_score = -obj_val / sum(demand_matrix,'all');
    % plotDemandFlow(network, transfer_graph_digraph, flow_solution);
end

%% Helper Functions

function transfer_graph = constructTransferGraph(network, event_set, max_changeover_time, train_capacity)
    %% Construct a graph of possible passenger/freight movements using train arrivals/departures
    n_stations = size(network.adjacency_matrix,1);
    n_stops = size(event_set, 2);

    % Initial Demand Source Nodes [1 : n_stations]
    % Train route nodes [n_stations + 1 : n_stations+n_stops]
    % Demand sink nodes [n_stations + n_stops + 1 : 2 * n_stations + n_stops]
    %         | sources| routes|  sinks|
    % sources |   0    |   A   |   I   |
    % routes  |   0    |   B   |   C   |
    % sinks   |   0    |   0   |   0   |
    transfer_graph = zeros(2 * n_stations + n_stops);

    for i_stop = 1:n_stops 
        % Connect demand source (A)
        transfer_graph(event_set(1, i_stop), n_stations + i_stop) = Inf;

        % Connect demand sink (C)
        transfer_graph(n_stations + i_stop, n_stations + n_stops + event_set(1,i_stop)) = Inf;

        % Add train trips as edges (B)
        if i_stop > 1
            if event_set(3, i_stop - 1) == event_set(3, i_stop) % same train
                transfer_graph(n_stations + i_stop - 1, n_stations + i_stop) = train_capacity;
            end
        end
    end

    % Staying at station (I)
    %transfer_graph(1:n_stations, n_stations + n_stops + 1:end) = eye(n_stations, n_stations) * Inf;
    %transfer_graph(isnan(transfer_graph)) = 0;

    for i_station = 1:n_stations
        % Add station changeovers as edges (B)
        % Find successive train visits within a changeover timeframe
        idx_station_visits = find(event_set(1,:) == i_station);
        station_visits = event_set(:, idx_station_visits);

        [~, idx_sort_station_idxs] = sort(station_visits(2,:));
        idx_sorted_station_visits = idx_station_visits(idx_sort_station_idxs);

        for i_idx_sorted_station_visits = 2:length(idx_sorted_station_visits)
            if event_set(2, idx_sorted_station_visits(i_idx_sorted_station_visits)) < event_set(2, idx_sorted_station_visits(i_idx_sorted_station_visits-1)) + max_changeover_time
                transfer_graph(n_stations + idx_sorted_station_visits(i_idx_sorted_station_visits-1), n_stations + idx_sorted_station_visits(i_idx_sorted_station_visits)) = Inf;
            end
        end
    end 
    % Connect demand sink (C)
    transfer_graph(n_stations + i_stop, n_stations + n_stops + event_set(1,i_stop)) = Inf;
end

function plotDemandFlow(network, transfer_graph_digraph, flow_solution)
    n_nodes = transfer_graph_digraph.numnodes;
    n_edges = transfer_graph_digraph.numedges;
    n_stations = size(network.adjacency_matrix,1);
    n_demands = n_stations^2 - n_stations;

    edge_traffic = zeros(n_edges, 1);
    for i_edge = 1:n_edges
        edge_traffic(i_edge) = sum(flow_solution(i_edge:n_edges:end-n_demands));
    end

    figure();
    heatmap = hot;
    colormap(heatmap(1:end-80,:));
    plot(transfer_graph_digraph, 'Layout', 'layered', 'Sources', [1:n_stations], 'Sinks', [n_nodes-n_stations+1:n_nodes], 'EdgeCData', edge_traffic, 'LineWidth', 2.5, 'MarkerSize', 5);
end

function distance = trainDistance(network, traj_set, i_train, j_train, timestep)
    %% Calculates distance of two trains given a set of trajectories
    % trajectory_set dimensions (n_trains, 3, timestep)
    % trajectory_set values (edge 0-n, position on edge 0-1, train orientation on edge -1,1)

    edge_i = int32(traj_set(i_train, 1, timestep));
    edge_j = int32(traj_set(j_train, 1, timestep));

    i_edge_length = network.edge_values(edge_i);
    j_edge_length = network.edge_values(edge_j);

    if edge_i == edge_j
        distance = abs(traj_set(i_train, 2, timestep) - traj_set(j_train, 2, timestep)) * i_edge_length;
    else
        % Find shortest path with precomputed distance matrix
        i_node_backward = network.edge_rows(edge_i);
        i_node_forward = network.edge_cols(edge_i);
        j_node_backward = network.edge_rows(edge_j);
        j_node_forward = network.edge_cols(edge_j);

        i_remaining_backward_length = traj_set(i_train, 2, timestep) * i_edge_length;
        i_remaining_forward_length = i_edge_length - i_remaining_backward_length;
        j_remaining_backward_length = traj_set(j_train, 2, timestep) * j_edge_length;
        j_remaining_forward_length = j_edge_length - j_remaining_backward_length;

        dist1 = i_remaining_backward_length + network.all_shortest_paths(i_node_backward,j_node_backward) + j_remaining_backward_length;
        dist2 = i_remaining_forward_length + network.all_shortest_paths(i_node_forward,j_node_backward) + j_remaining_backward_length;
        dist3 = i_remaining_backward_length + network.all_shortest_paths(i_node_backward,j_node_forward) + j_remaining_forward_length;
        dist4 = i_remaining_forward_length + network.all_shortest_paths(i_node_forward,j_node_forward) + j_remaining_forward_length;

        distance = min([dist1, dist2, dist3, dist4]);
    end
end

%% Solution Generation 

function solution = randomSolution(params)
    solution = rand(params.n_trains, params.n_timesteps * 2);
end

function solution = greedySolution(network, params)
    %% Brute-forces valid routes sequentially
    solution = rand(1, params.n_timesteps * 2);
    traj_set(1,:,:) = constructTrajectory(network, solution, params.initial_positions(1,:), params.initial_speeds(1), params.max_accel, params.max_speed);
    i_train = 2;
    while i_train <= params.n_trains
        disp(["Greedy Placement: ", mat2str(round(i_train/params.n_trains * 100)), "%"]);
        abort = 1;
        collision_score = -Inf;
        while collision_score < 0
            new_train_solution = rand(1, params.n_timesteps * 2);
            new_traj = constructTrajectory(network, new_train_solution, params.initial_positions(i_train,:), params.initial_speeds(i_train), params.max_accel, params.max_speed);
            new_traj_set = cat(1, traj_set, reshape(new_traj, 1, size(new_traj,1), size(new_traj,2)));
            new_solution = cat(1, solution, new_train_solution);
            collision_score = collisionPenalties(network, new_traj_set, params.min_separation, params.max_speed);

            abort = abort + 1;
            if abort > 1000
                warning("Failed to find valid placement. Trying from scratch..");
                new_solution = rand(1, params.n_timesteps * 2);
                new_traj_set = [];
                new_traj_set(1,:,:) = constructTrajectory(network, new_solution, params.initial_positions(1,:), params.initial_speeds(1), params.max_accel, params.max_speed);
                i_train = 1;
                break;
            end
        end
        traj_set = new_traj_set;
        solution = new_solution;
        i_train = i_train + 1;
    end
end

%% Search Methods

function greedyRandomSearch(network, params)
    %% Generate a random collision free solution then evaluate demand score
    csvwrite("network.csv", network.adjacency_matrix);
    csvwrite("demands.csv", params.demand_matrix);
    tic;
    best_solution_set = {0 -Inf 0 0}; % traj_set, demand_score, transfer_graph_digraph, flow_solution
    while best_solution_set{2} < 1
        [traj_set, event_set] = constructTrajectorySet(network, greedySolution(network, params), params.initial_positions, params.initial_speeds, params.max_accel, params.max_speed);
        [new_demand_score, new_transfer_graph_digraph, new_flow_solution] = demandSatisfaction(network, event_set, params.demand_matrix, params.max_changeover_time, params.train_capacity);
        if best_solution_set{2} < new_demand_score
            best_solution_set = {traj_set new_demand_score new_transfer_graph_digraph new_flow_solution}
            csvwrite("trajectories_edges.csv", squeeze(traj_set(:,1,:)));
            csvwrite("trajectories_positions.csv", squeeze(traj_set(:,2,:)));
        end
        if toc > 30
            break;
        end
    end
    plotDemandFlow(network, best_solution_set{3}, best_solution_set{4});
end

greedyRandomSearch(network, params);