function transfer_graph = constructTransferGraph(network, event_set, max_changeover_time, train_capacity)
    %% Construct a graph of possible passenger/freight movements using train arrivals/departures
    n_stations = size(network.adjacency_matrix,1);
    n_stops = size(event_set, 1);

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
        transfer_graph(event_set(i_stop, 3), n_stations + i_stop) = Inf;

        % Connect demand sink (C)
        transfer_graph(n_stations + i_stop, n_stations + n_stops + event_set(i_stop, 3)) = Inf;

        % Add train trips as edges (B)
        if i_stop > 1
            if event_set(i_stop - 1, 1) == event_set(i_stop, 1) % same train
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
        idx_station_visits = find(event_set(:, 3) == i_station);
        station_visits = event_set(idx_station_visits, :);

        [~, idx_sort_station_idxs] = sort(station_visits(2,:));
        idx_sorted_station_visits = idx_station_visits(idx_sort_station_idxs);

        for i_idx_sorted_station_visits = 2:length(idx_sorted_station_visits)
            if event_set( idx_sorted_station_visits(i_idx_sorted_station_visits), 2) < event_set(idx_sorted_station_visits(i_idx_sorted_station_visits-1), 2) + max_changeover_time
                transfer_graph(n_stations + idx_sorted_station_visits(i_idx_sorted_station_visits-1), n_stations + idx_sorted_station_visits(i_idx_sorted_station_visits)) = Inf;
            end
        end
    end 
    % Connect demand sink (C)
    transfer_graph(n_stations + i_stop, n_stations + n_stops + event_set(i_stop, 3)) = Inf;
end

