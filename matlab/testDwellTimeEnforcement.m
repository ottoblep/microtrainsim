timesteps = 100;
interp_steps = 10;
points = rand(1, interp_steps) * timesteps;
vals = rand(1, interp_steps);
a_max = 0.1;
v_max = 1;
v_init = 2 * (rand(1,1) - 0.5) * v_max;
dwell_time = 0.1 * timesteps;

acceleration = (interpolateSolutionCurve(points, vals, 1:timesteps) * 2 - 1) * a_max;
speeds = v_init + cumtrapz(acceleration);
acceleration(speeds>v_max & acceleration>0) = 0;
acceleration(speeds<-v_max & acceleration<0) = 0;
speeds = v_init + cumtrapz(acceleration);
position = cumtrapz(speeds);

initial_arrival_time = 0.5 * timesteps;

% Braking at max acceleration until standstil leads to a quadratic position curve
% We calculate the timestep to start braking so the train comes to a halt exactly at the stop 
% Formula for the start of braking
% speed = a_max * sqrt(2*d / a_max)
% d distance to stop, speed current_speed, a_max acceleration

stop_position = position(initial_arrival_time);
approach_direction = sign(speeds(initial_arrival_time));

pre_stop_timesteps = 1:initial_arrival_time;
first_approach_idx = find(sign(speeds(1:initial_arrival_time)) ~= approach_direction, 1, 'last') + 1;
if isempty(first_approach_idx)
    first_approach_idx = 1;
end
approach_timesteps = first_approach_idx:initial_arrival_time;

% Determine ideal start of braking
d = abs(position(approach_timesteps) - stop_position);
start_braking_timestep = first_approach_idx + find(d>=0.5 * speeds(approach_timesteps).^2 / a_max, 1, 'last');
% Pessimistic rounding for discretization error
new_arrival_time = ceil(start_braking_timestep + (abs(speeds(start_braking_timestep) / a_max)));
exact_acceleration = abs(speeds(start_braking_timestep) / (new_arrival_time - start_braking_timestep + 1));

% Modify acceleration curve for ideal approach
acceleration(start_braking_timestep: new_arrival_time) = -approach_direction * exact_acceleration;
acceleration(new_arrival_time:new_arrival_time + dwell_time) = 0;
speeds = v_init + cumtrapz(acceleration);
acceleration(speeds>v_max & acceleration>0) = 0;
acceleration(speeds<-v_max & acceleration<0) = 0;
speeds = v_init + cumtrapz(acceleration);
%speeds(new_arrival_time:new_arrival_time + dwell_time) = 0;
position = cumtrapz(speeds);

v_error = speeds(new_arrival_time)
p_error = position(new_arrival_time) - stop_position
startpoint_error = min(abs(0.5 * speeds(approach_timesteps).^2 / a_max - d))

clf;
hold on;
%scatter(points, vals);
scatter(first_approach_idx, 0,'DisplayName','first\_approach\_idx');
scatter(initial_arrival_time, 0,'DisplayName','initial\_arrival\_time');
scatter(start_braking_timestep, 0,'DisplayName','start\_braking\_timestep');
scatter(new_arrival_time, 0,'DisplayName','new\_arrival\_time');
plot(1:timesteps, acceleration,'DisplayName','accel');
plot(1:timesteps, speeds,'DisplayName','speed');
plot(1:timesteps, position - stop_position,'DisplayName','dist\_to\_stop');
legend();

function y_new = interpolateSolutionCurve(x, y, x_new)
    %% Interpolate sparse curve representation to continuous one and normalize
    [~ , unique_idxs, ~] = unique(x);
    y_new = interp1(x(unique_idxs), y(unique_idxs), x_new, 'linear', 'extrap');
    y_new(y_new>1) = 1;
    y_new(y_new<0) = 0;
end