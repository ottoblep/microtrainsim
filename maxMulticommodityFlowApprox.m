%% Max Multicommodity Network Flow algorithm from "Approximating fractional multicommodity flow independent of the number of commodities" (Fleischer 2000)
% This implementation is adapted to a directed graph without loops and and and identical number of source and sink nodes
function [flow_value, edge_flows] = maxMulticommodityFlowApprox(network, network_digraph, n_source_sink_nodes, demand_matrix, e_accuracy)
     % network is the directed adjacency matrix with format (sources, other, sinks)
     % e_accuracy is the epsilon in (1+epsilon)

     k_demand_pairs = n_source_sink_nodes^2 - n_source_sink_nodes;

     % We add helper nodes to conform to limit flow for each source-sink-pair according to the demand matrix
     % This adds n_source_sink_nodes^2 - n_source_sink_nodes nodes to the network with each one edge that has capacity demand_i,j

     % Example for three source and sink nodes

     %         | demand_helper_sources       | sources| routes|  sinks|
     % helpers |                             |D12     |       |       |
     %         |                             |D13     |       |       |
     %         |             0               |   D21  |       |       |
     %         |                             |   D23  |   0   |   0   |
     %         |                             |     D31|       |       |
     %         |                             |     D32|       |       |
     % sources |                             |   0    |   A   |   I   |
     % routes  |             0               |   0    |   B   |   C   |
     % sinks   |                             |   0    |   0   |   0   |

     [demand_rows,demand_cols,demand_vals] = find(demand_matrix);
     k_demand_pairs = numel(demand_vals);
     n_helper_nodes = k_demand_pairs;

     % Insert transfer graph (A,B,C,I)  
     network_full = sparse(size(network,1) + n_helper_nodes, size(network,1) + n_helper_nodes);
     network_full(n_helper_nodes + 1:n_helper_nodes + size(network,1), n_helper_nodes + 1:n_helper_nodes + size(network,1)) = network;

     % Add demand constraint helper nodes (D12, D13 etc)
     for i_demand_pair = 1:n_helper_nodes
          network_full(i_demand_pair, n_helper_nodes + demand_rows(i_demand_pair)) = demand_vals(i_demand_pair);
     end

     network_full_digraph = digraph(network_full);
     n_nodes = size(network_full,1);
     m_edges = numedges(network_full_digraph);
     edge_capacities = network_full_digraph.Edges.Weight;

     demand_paths = cell(k_demand_pairs, 1);
     for i_demand_pair = 1:k_demand_pairs
          [~, demand_paths{i_demand_pair}] = allpaths(network_full_digraph, i_demand_pair, n_nodes - n_source_sink_nodes + demand_cols(i_demand_pair));
     end
     n_demand_path_sizes = cellfun('size', demand_paths, 1);
     n_demand_paths = sum(n_demand_path_sizes);

     % Fully Polynomial-Time Approximation Scheme Algorithm

     % L maximum number of arcs in augmenting path
     L = max(n_demand_path_sizes);

     % delta is initial dual problem path length
     delta = (1 + e_accuracy) / ((1 + e_accuracy) * L)^(1/e_accuracy);

     % Initialize decision vars for primal and dual problem
     x = zeros(1, n_demand_paths);
     l = ones(1, m_edges) * delta;

     n_phases = 0;
     for i = 1:floor(log(((1+e_accuracy)/delta)) / log(1+e_accuracy))
          n_phases = n_phases + 1;
          for j = 1:k_demand_pairs
               relevant_paths = demand_paths{j};
               if isempty(relevant_paths)
                    continue;
               end
               [P_idx, P_len] = weightedShortestPath(relevant_paths, l);
               P = relevant_paths{P_idx,1};
               global_path_idx_base = sum(n_demand_path_sizes(1:j-1));

               while P_len < min([1 delta*(1+e_accuracy)^i])
                    demand_to_assign_u = min(edge_capacities(P));

                    global_path_idx = P_idx + global_path_idx_base;
                    x(global_path_idx) = x(global_path_idx) + demand_to_assign_u;

                    l(P) = l(P) + ((l(P) * demand_to_assign_u * e_accuracy) ./ edge_capacities(P)');

                    [P_idx, P_len] = weightedShortestPath(relevant_paths, l);
                    P = relevant_paths{P_idx,1};
               end
          end
     end

     % Scale result
     result = x / (log(1 / delta) / log(1 + e_accuracy));
     flow_value = sum(result);

     % Sum edge flows
     edge_flows = zeros(m_edges, 1);
     for j_demand_pair = 1:k_demand_pairs
          paths = demand_paths{j_demand_pair};
          if isempty(paths)
               continue
          end
          for i_path = 1:numel(paths)
               P = paths{i_path,1};
               global_path_idx = i_path + sum(n_demand_path_sizes(1:j_demand_pair-1));
               for i_path_edge = 1:numel(P)
                    global_edge_idx = P(i_path_edge);
                    edge_flows(global_edge_idx) = edge_flows(global_edge_idx) + result(global_path_idx);
               end
          end
     end

     % Cut helper edges
     edge_flows = edge_flows(end - numedges(network_digraph) + 1:end);
end

function [shortest_path_idx best_length] = weightedShortestPath(paths, edge_weights)
     best_length = Inf;
     shortest_path_idx = 0;
     for i_path = 1:numel(paths)
          path = paths{i_path,1};
          length = sum(edge_weights(path));

          if length < best_length
               best_length = length;
               shortest_path_idx = i_path;
          end
     end
end