transfer_graph = [0 0 0 5 0 0 0 0;
                  0 0 0 5 5 0 0 0;
                  0 0 0 0 5 0 0 0;
                  0 0 0 0 5 5 5 0;
                  0 0 0 0 0 0 5 5;
                  0 0 0 0 0 0 0 0;
                  0 0 0 0 0 0 0 0;
                  0 0 0 0 0 0 0 0;];
transfer_graph_digraph = digraph(transfer_graph);

demand_matrix = [0 3 0;
                 0 0 0;
                 0 0 0;];

e_accuracy = 0.1;

tic
[flow_value1, edge_flows(1,:)] = maxMulticommodityFlowApprox(transfer_graph, transfer_graph_digraph, size(demand_matrix,1), demand_matrix, e_accuracy);
toc
tic
[flow_value2, edge_flows(2,:)] = maxMulticommodityFlowLP(transfer_graph, transfer_graph_digraph, size(demand_matrix,1), demand_matrix);
toc
%plotDemandFlow(transfer_graph_digraph, size(demand_matrix,1), edge_flows(1,:));