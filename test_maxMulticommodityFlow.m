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
                  0 15 0;];

e_accuracy = 0.01;

n_source_sink_nodes = 3;

[flow_value, edge_flows] = maxMulticommodityFlowApprox(transfer_graph, transfer_graph_digraph, size(demand_matrix,1), demand_matrix, e_accuracy);
[flow_value, edge_flows] = maxMulticommodityFlowLP(transfer_graph, transfer_graph_digraph, size(demand_matrix,1), demand_matrix);
