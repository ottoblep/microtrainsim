timesteps = 150;
interp_steps = 5;
points = rand(1, interp_steps) * timesteps;
vals = rand(1, interp_steps);
a_max = 0.1;
v_max =  1;
x_0 = rand(1,1);
v_init = 2 * (rand(1,1) - 0.5) * v_max;
dwell_time = round(0.3 * timesteps);

close all;
hold on;

acceleration = (interpolateSolutionCurve(points, vals, 1:timesteps) * 2 - 1) * a_max;
speeds = v_init + cumsum(acceleration);
speeds(speeds>v_max) = v_max;
speeds(speeds<-v_max) = -v_max;
position = x_0 + cumsum(speeds);

initial_arrival_time = round(0.5 * timesteps);

approach_direction = sign(speeds(initial_arrival_time));

pre_stop_timesteps = 1:initial_arrival_time;
first_approach_idx = find(sign(speeds(1:initial_arrival_time)) ~= approach_direction, 1, 'last') + 1;
if isempty(first_approach_idx)
    first_approach_idx = 1;
end
approach_timesteps = first_approach_idx:initial_arrival_time;

% Discrete Formula for distance covered under k braking steps
% k*v(n) - a*k*(k+1)/2 
% maximum (stop point) at (v-a/2)/a

% Select timestep to start braking (will always overshoot slightly)
v = speeds(approach_timesteps);
k = floor(v./a_max);
d_brake = k .* v - a_max * k .* (k+1) * 0.5;
d = abs(position(approach_timesteps) - position(initial_arrival_time));
start_braking_timestep = first_approach_idx + find(d < d_brake, 1, 'first') + 1;

% Calculate braking curve
n_full_braking_steps = floor(abs(speeds(start_braking_timestep-1))/a_max);
if n_full_braking_steps + 1 > timesteps - start_braking_timestep
    error("Not enough time to brake.");
end
a_last_step = abs(speeds(start_braking_timestep-1)) - n_full_braking_steps * a_max;
assert((n_full_braking_steps * a_max + a_last_step == abs(speeds(start_braking_timestep-1))));
end_braking_timestep = start_braking_timestep + n_full_braking_steps;

scatter(initial_arrival_time, position(initial_arrival_time),'DisplayName','initial arrival time');
plot(1:timesteps, position,'DisplayName','position old');

% Adjust acceleration values
acceleration = [speeds(1)-v_init diff(speeds)];
acceleration(start_braking_timestep:end_braking_timestep - 1) = -approach_direction * a_max;
acceleration(end_braking_timestep) = -approach_direction * a_last_step;
if end_braking_timestep + dwell_time > timesteps
    acceleration(end_braking_timestep+1:timesteps) = 0;
else
    acceleration(end_braking_timestep+1:end_braking_timestep + dwell_time) = 0;
end
% Recalculate positions
speeds = v_init + cumsum(acceleration);
position = x_0 + cumsum(speeds);

v_error = speeds(end_braking_timestep)
p_relative_error = approach_direction * (position(end_braking_timestep) - position(initial_arrival_time))

scatter(first_approach_idx, position(first_approach_idx),'DisplayName','first approach idx');
scatter(start_braking_timestep, position(start_braking_timestep), 'DisplayName','start braking timestep');
scatter(end_braking_timestep, position(end_braking_timestep),'DisplayName','end braking timestep');
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