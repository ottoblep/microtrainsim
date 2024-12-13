function plotSpeedEnvelope(network, params, v_targets, traj)
    clf; hold on;
    plot(traj(:,4));
    plot(network.speed_limits(traj(:,1)));
    plot(-network.speed_limits(traj(:,1)));
    scatter(v_targets(:, 1), v_targets(:, 2));
end