%% Link-Flow Model Linear Programming Formulation of the Multicommodity Flow Problem
% This implementation is adapted to a directed graph without loops and and and identical number of source and sink nodes
function [flow_value, edge_flows] = maxMulticommodityFlowLP(network, network_digraph, n_stations, demand_matrix)
    % network is the directed adjacency matrix

    network_edge_idxs = find(network);
    n_nodes = size(network,1); % n_stations + n_stops + n_stations
    n_edges = numel(network_edge_idxs);
    % Only nondiagonal entries in the demand matrix
    n_demands = n_stations^2 - n_stations;
    n_decision_vars = (n_edges + 1)* n_demands;

    %% Splittable Multi-commodity maximum flow problem https://en.wikipedia.org/wiki/Multi-commodity_flow_problem
    % max MCNF maximizes global flow
    % max-concurrent MCNF maximizes the satisfied percentage of all flows
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
        out_edges_idxs = outedges(network_digraph, i_node);
        in_edges_idxs = inedges(network_digraph, i_node);
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
    b(1:n_edges*n_demands) = repmat(network(network_edge_idxs), n_demands, 1);
    demands_transposed = transpose(demand_matrix);
    b(end - n_demands + 1:end) = demands_transposed(~eye(size(demand_matrix)));
    lb = sparse(n_decision_vars, 1);
    options = optimoptions('linprog','Display','none');
    % Run Solver
    [flow_solution, obj_val, ~, output] = linprog(f, A, b, Aeq, beq, lb, [], options);

    edge_flows = zeros(n_edges, 1);
    for i_edge = 1:n_edges
        edge_flows(i_edge) = sum(flow_solution(i_edge:n_edges:end-n_demands));
    end

    flow_value = -obj_val;
end