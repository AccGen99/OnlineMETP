function val = run(tf)
clc
Params_Calcs;
vals = merv_ode_calcs(tf, tau_w);
% Lagragian for x = [x, y, phi] and x_dot = [x_dot, y_dot, phi_dot]
alpha = [0, 0, 0]; %[alpha_x, alpha_y, alpha_phi];   x
lambda = [0, 0, 0]; %[lambda_xdot, lambda_ydot, lambda_phidot];  x_dot

% Vector of phi_dot values
t_span = linspace(0, tf, 100);
phi_dot = 0*t_span;  % Array vector, phi_dot(i) corresponds to time(i)
range = 1:length(t_span);
h2 = vals(1);
h1 = vals(2);
% Rotational Velocity
for i=range
    phi_dot(i) = -MERV(phi_f, t_span(i), tf, tau_w, h2, h1)/(2*pi);
end
%phi_final_calc = (t_span(2)-t_span(1))*trapz(phi_dot)
figure
plot(t_span, phi_dot)
title('Angular Velocity vs Time')
xlabel('Time')
ylabel('Angular Velocity')
% alpha_phi
alpha(3) = (-1*(C3^1.5)*phi_f/(C4*(C3^0.5)*tf - 2*tanh((C3^0.5)*tf/2)))/(2*pi);

% Translational velocity in x-direction
% Requires appropriate values of initial acceleration and alpha_x
% For positive xf, alpha_x < 0 and a_x(0) > 0
u_max = [0, -1, 1]';
ax_max = maxAcc(k2, R, A, B, u_max); % Max possible initial acceleration

%alpha_x = -4; % Starting with an arbitrary negative number
tspan = t_span;%[0, tf];
ax_0 = 0; % Initial calculation for zero initial-acceleration condition
ic = [0, ax_0]; % Ini_vel & ini_acc
opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
[t_1,vx_1]=ode45(@(t,y) vx_ode(t, y, t_span, phi_dot, C1, C2, alpha_x),...
    tspan,ic,opts);
vx_tf = vx_1(end,1);

% From linearity condition
fL = vx_tf;

ic2 = [0, ax_max];
[t_2,vx_2] = ode45(@(t,y) vx_ode(t, y, t_span, phi_dot, C1, C2,...
    alpha_x), tspan, ic2, opts);
vx_tf2 = vx_2(end,1);

% From linearity condition
fU = vx_tf2;
ax_0_t = ax_max*(fL/(fL-fU));

% Finding final position candidate for x coordinate
ic3 = [0, ax_0_t];
[t_3,vx_3] = ode45(@(t,y) vx_ode(t, y, t_span, phi_dot, C1, C2,...
    alpha_x), tspan, ic3, opts);
t_step = t_3(2) - t_3(1);
vx_vector = vx_3(:,1);
ax_vector = vx_3(:,2);
% Integrating velocity & acceleration assuming constant acceleration in
% each time interval
rng = 1:length(t_3);
c_xf = 0;
for j = rng
    ini_vel = vx_vector(j);
    accln = ax_vector(j);
    c_xf = c_xf + ini_vel*t_step + 0.5*accln*t_step^2;
end
c_xf;
s_factor = c_xf/xf;%disp/xf;
alphax_scaled = alpha_x*C2/(s_factor);%C2*s_factor*alpha_x
ax_0_t_scaled = ax_0_t/s_factor;%s_factor*ax_0_t

alpha(1) = alphax_scaled;

% Validating
ic4 = [0, ax_0_t_scaled];
[t_4,vx_4] = ode45(@(t,y) vx_ode(t, y, t_span, phi_dot, C1, C2,...
    alphax_scaled), tspan, ic4, opts);
t_step = t_4(2) - t_4(1);
vx_vector1 = vx_4(:,1);
ax_vector1 = vx_4(:,2);

figure
plot(t_span, vx_4(:,1))
title('X Velocity vs Time')
xlabel('Time')
ylabel('X Velocity')


rng = 1:length(t_4);
disp = 0;
for j = rng
    ini_vel = vx_vector1(j);
    accln = ax_vector1(j);
    disp = disp + ini_vel*t_step + 0.5*accln*t_step^2;
end
disp
% Velocity for x, y and phi
vx_vectOr = interp1(t_4, vx_vector1, t_span);
vy_vector = 0*vx_vectOr;
ang_vel = phi_dot;

% Acceleration for x, y and phi
ax_vectOr = interp1(t_4, ax_vector1, t_span);
ay_vector = 0*ax_vectOr;
ang_acc = 0*phi_dot;
range1 = 1:length(phi_dot);
for k = range1
    ang_acc(k) = C3*phi_dot(k) + C4*alpha(3);
end
x_dot_dot = [ax_vectOr', ay_vector', ang_acc'];
x_dot = [vx_vectOr', vy_vector', ang_vel'];

% Finding lambda
lamda = 0*x_dot;
for l = range1
    xdot_vec = x_dot(l,:);
    xdotdot_vec = x_dot_dot(l,:);
    omega = ang_vel(l);
    R_dot = omega*R_dot;
    lamda(l,:) = FindLambda(w2, A, xdot_vec', k2, w1, Q, xdotdot_vec',...
        C, omega, R, R_dot);
end

% Finding lamda dot
lamda_dot = DerLamda(lamda, tf);

% Find heading angle
u = 0*lamda_dot;
t_step0 = t_span(2) - t_span(1);
phi = phiCalc(phi_dot, ang_acc, t_step0);

% Find control vector
for m = range1
    l_dot = lamda_dot(m,:);
    omega = ang_vel(m);
    p = phi(m);
    R_d = R_dot_calc(p);
    l = lamda(m,:);
    u(m,:) = cVec(B, R, l_dot', alpha', omega, R_d, l', C, A, w2);
end

E = 0*range1;
for n=range1
    u_ = u(n,:);
    x_dot_=x_dot(n,:);
    p = phi(n);
    R = R_calc(p);
    E(n) = energyCalc(w1,u_',w2,x_dot_',R,B);
end
Energy = t_step*trapz(E);
Energy
u_max = 1;
u_max_obtd = max(u);
u_max_obtd = max(u_max_obtd);
u_max_obtd = abs(u_max_obtd);

u_min_obtd = min(u);
u_min_obtd = min(u_min_obtd);
u_min_obtd = abs(u_min_obtd);

if u_min_obtd > u_max_obtd
    u_max_obtd = u_min_obtd;
end

val = u_max_obtd - u_max;
end
