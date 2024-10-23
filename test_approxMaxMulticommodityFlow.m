network = [0 0 0 5 0 0 0 0;
           0 0 0 6 5 0 0 0;
           0 0 0 0 8 0 0 0;
           0 0 0 0 5 5 5 0;
           0 0 0 0 0 0 5 5;
           0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0;
           0 0 0 0 0 0 0 0;];

demands = [10 0 7 0 0 15];

e_accuracy = 0.01;

n_source_sink_nodes = 3;

result = approxMaxMulticommodityFlow(network, n_source_sink_nodes, demands, e_accuracy);