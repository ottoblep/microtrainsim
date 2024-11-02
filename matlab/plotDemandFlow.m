function plotDemandFlow(transfer_graph_digraph, n_source_sink_nodes, edge_flows)
    n_nodes = transfer_graph_digraph.numnodes;
    n_edges = transfer_graph_digraph.numedges;
    n_stations = n_source_sink_nodes;
    n_demands = n_stations^2 - n_stations;

    figure();
    heatmap = hot;
    colormap(heatmap(1:end-80,:));
    plot(transfer_graph_digraph, 'Layout', 'layered', 'Sources', [1:n_stations], 'Sinks', [n_nodes-n_stations+1:n_nodes], 'EdgeCData', edge_flows, 'LineWidth', 2.5, 'MarkerSize', 5);
end