transfer_graph = [0 0 0 5 0 0 0 0;
                  0 0 0 6 5 0 0 0;
                  0 0 0 0 8 0 0 0;
                  0 0 0 0 5 5 5 0;
                  0 0 0 0 0 0 5 5;
                  0 0 0 0 0 0 0 0;
                  0 0 0 0 0 0 0 0;
                  0 0 0 0 0 0 0 0;];
transfer_graph_digraph = digraph(transfer_graph);

demand_matrix = [0 10 0;
                  7  0 0;
                  3 15 0;];

e_accuracy = 0.01;

[flow_value1, edge_flows1] = maxMulticommodityFlowApprox(transfer_graph, transfer_graph_digraph, size(demand_matrix,1), demand_matrix, e_accuracy);%
[flow_value2, edge_flows2] = maxMulticommodityFlowLP(transfer_graph, transfer_graph_digraph, size(demand_matrix,1), demand_matrix);