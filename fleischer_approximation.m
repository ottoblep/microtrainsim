%% Max-concurrent Multicommodity Network Flow algorithm from "Approximating fractional multicommodity flow independent of the number of commodities" (Fleischer 2000)

% Only relevant demands are listed (no diagonal and only upper) e.g. 1-7, 1-8, 2-6, 2-8, 3-6, 3-7,
demands = [10 0 7 0 0 15];

network = [0 0 0 5 0 0 0 0;
           0 0 0 6 5 0 0 0;
           0 0 0 0 8 0 0 0;
           0 0 0 0 5 5 5 0;
           0 0 0 0 0 0 5 5;
           0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0;];

g_network = digraph(network);
[~, ~, edge_capacities] = find(network);

n_nodes = size(network,1);
m_edges = numel(edge_capacities);
n_flows = 3;
k_demand_pairs = n_flows^2 - n_flows;
demand_pair_idxs = combinations(1:n_flows, n_nodes-n_flows+1:n_nodes);
% Remove diagonal
demand_pair_idxs = demand_pair_idxs{demand_pair_idxs{:,1}~=demand_pair_idxs{:,2} - (n_nodes-n_flows),:};

demand_paths = cell(k_demand_pairs, 1);
for i_demand_pair = 1:k_demand_pairs
     demand_paths{i_demand_pair} = allpaths(g_network, demand_pair_idxs(i_demand_pair, 1), demand_pair_idxs(i_demand_pair, 2));
end
n_demand_path_sizes = cellfun('size', demand_paths, 2);
n_demand_paths = sum(n_demand_path_sizes);
e_accuracy = 0.01;

% delta is initial dual problem path length
delta = (m_edges / (1 - e_accuracy))^(-1/e_accuracy);

% Initialize decision vars for primal and dual problem
x = zeros(1, n_demand_paths);
l = delta / edge_capacities;

% Dual objective function D(l) = sum(u(e)l(e))
D_l = @(l_arg) sum(edge_capacities.*l_arg');
demand_remaining = zeros(1,k_demand_pairs);

n_phases = 0;
while D_l(l) < 1
     n_phases = n_phases + 1;
     for j = 1:k_demand_pairs
          demand_remaining(j) = demands(j);
          while D_l(l) < 1 & demand_remaining(j) > 0
               % Shortest path in P_j using l
               relevant_paths = paths{j};
               P_idx = weightedShortestPath(relevant_paths, l);
               P = relevant_paths{1,P_idx};

               smallest_capacity_in_p = min(edge_capacities(P));
               demand_to_assign_u = min([demand_remaining(j) smallest_capacity_in_p]);

               demand_remaining(j) = demand_remaining(j) - demand_to_assign_u;

               global_path_idx = P_idx + sum(n_demand_path_sizes(1:j-1));

               % Assign demand
               x(global_path_idx) = x(global_path_idx) + demand_to_assign_u;

               for i_path_edge = 1:numel(P)
                    global_edge_idx = P(i_path_edge);
                    l(global_edge_idx) = l(global_edge_idx) * (1 + (demand_to_assign_u * e_accuracy) / edge_capacities(global_edge_idx));
               end
          end
     end
end

% Scale result
result = x / (log(1 / delta) / log(1 + e_accuracy));

function shortest_path_idx = weightedShortestPath(paths, edge_weights)
     best_length = Inf;
     for i_path = 1:numel(paths)
          path = paths{1,i_path};
          length = 0;
          for j_edge_in_path = 1:size(path,2)
               length = length + edge_weights(path(j_edge_in_path));
          end
          if length < best_length
               best_length = length;
               shortest_path_idx = i_path;
          end
     end
end