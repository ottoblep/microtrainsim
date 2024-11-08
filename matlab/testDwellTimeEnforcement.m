timesteps = 20;
interp_steps = 3;
points = rand(1, interp_steps) * timesteps;
vals = rand(1, interp_steps);
a_max = 0.1;
v_max = 1;
x_0 = rand(1,1);
v_init = 2 * (rand(1,1) - 0.5) * v_max;
dwell_time = round(0.2 * timesteps);

clf;
hold on;

acceleration = (interpolateSolutionCurve(points, vals, 1:timesteps) * 2 - 1) * a_max;
speeds = v_init + cumtrapz(acceleration);
acceleration(speeds>v_max & acceleration>0) = 0;
acceleration(speeds<-v_max & acceleration<0) = 0;
speeds = v_init + cumtrapz(acceleration);
position = x_0 + cumtrapz(speeds);
plot(1:timesteps, position,'DisplayName','position\_old');

initial_arrival_time = round(0.5 * timesteps);
scatter(initial_arrival_time, position(initial_arrival_time),'DisplayName','initial\_arrival\_time');

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
start_braking_timestep = first_approach_idx - 1 + find(d>=0.5 * speeds(approach_timesteps).^2 / a_max, 1, 'last');
scatter(start_braking_timestep, position(start_braking_timestep), 'DisplayName','start\_braking\_timestep');
% Pessimistic rounding for discretization error
new_arrival_time = start_braking_timestep + ceil(abs(speeds(start_braking_timestep) / a_max)) + 1;
exact_acceleration = abs(speeds(start_braking_timestep)) / (abs(new_arrival_time - start_braking_timestep - 1));

% Modify acceleration curve for ideal approach
acceleration(start_braking_timestep + 1:new_arrival_time - 1) = -approach_direction * exact_acceleration;
acceleration(new_arrival_time:new_arrival_time + dwell_time) = 0;
speeds = v_init + cumtrapz(acceleration);1
acceleration(speeds>v_max & acceleration>0) = 0;
acceleration(speeds<-v_max & acceleration<0) = 0;
speeds = v_init + cumtrapz(acceleration);
%speeds(new_arrival_time:new_arrival_time + dwell_time - 1) = 0;
%assert(all(abs(diff(speeds)) < a_max));
position = x_0 + cumtrapz(speeds);

v_error = speeds(new_arrival_time)
p_error = position(new_arrival_time) - stop_position
startpoint_error = min(abs(0.5 * speeds(approach_timesteps).^2 / a_max - d))

scatter(first_approach_idx, position(first_approach_idx),'DisplayName','first\_approach\_idx');

scatter(new_arrival_time, position(new_arrival_time),'DisplayName','new\_arrival\_time');
plot(1:timesteps, acceleration,'DisplayName','accel');
plot(1:timesteps, speeds,'DisplayName','speed');
plot(1:timesteps, position,'DisplayName','position\_new');
legend();

function y_new = interpolateSolutionCurve(x, y, x_new)
    %% Interpolate sparse curve representation to continuous one and normalize
    [~ , unique_idxs, ~] = unique(x);
    y_new = interp1(x(unique_idxs), y(unique_idxs), x_new, 'linear', 'extrap');
    y_new(y_new>1) = 1;
    y_new(y_new<0) = 0;
end