clear
clc
sysParams;

alpha_x = -1.0; % Arbitrary number

tspan = [0 tf];
y0 = [0 0];
[t,y_ini] = ode45(@(t,y) metv(t, y, alpha_x), tspan, y0); % Assuming initial velocity is 0

f_L = y_ini(end, 1);

max_a_x = accMax();

y0_0 = [0 max_a_x];
[t,y_max] = ode45(@(t,y) metv(t, y, alpha_x), tspan, y0_0); % Assuming initial velocity is 0

f_U = y_max(end, 1);

a_x_0_true = max_a_x*f_L/(f_L-f_U);
y_0_0 = [0 a_x_0_true];
[t,y_true] = ode45(@(t,y) metv(t, y, alpha_x), tspan, y_0_0); % Assuming initial velocity is 0

velocityArray = y_true(:,1);
c_x_f = sum(velocityArray)*t(end);
scaling_factor = c_x_f/x_f;

alpha_x = alpha_x*scaling_factor;
a_x_0_true = a_x_0_true*scaling_factor;