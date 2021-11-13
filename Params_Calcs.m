% Initial & Final state values
phi = [0, pi]; % phi(1) is initial and phi(2) is final
xf = 0.3; % Final x-coordinate
tf = 0.89655;%1.0;    % Apparently we're supposed to input this :/
zi = [0, 0, phi(1), 0, 0, 0]';
zf = [xf, 0, phi(2), 0, 0, 0]';

% System parameter values
Fv = 0.0025; % (N.m/(rad/sec))
Kb = 0.0255; % (V/(rad/sec))
Kt = 0.0255; % (N.m/A)
Mc = 1.59;      % Kg                Mass of robot body
Mw = 0.133;     % Kg                Mass of one wheel
L = 0.1255;     % m
Vs = 24;        % V
r = 0.03;       % m
tw = 0.045;     % m                 Nobody got no idea what it is/does
Ra = 7.9;       % ohm
n = 1;          % Gear ratio - motor & wheel

% Calculation of some constants
M = Mc + 3*Mw;
Jw = 0.5*Mw*r^2;
Jc = 0.5*Mc*L^2;
J = 3*Jw + Jc;
B = [0, 1, L; -sin(pi/3), -cos(pi/3), L; sin(pi/3), -cos(pi/3), L];
M = diag([M, M, J]);
k2 = Vs*Kt*n/(Ra*r);
D = k2*B';
k1 = 1.5*(Fv+Kt*Kb*n^2/Ra)/r^2;
C = k1*diag([1, 1, 2*L^2]);
Q = B'*B;
A = M+Jw*Q/r^2;
w1 = (Vs^2)/Ra;
w2 = Kb*n*Vs/(Ra*r);

% For differential equation of phi-dot
phi_f = phi(2);
C3 = (C(3,3)/A(3,3))^2 - k2*w2*Q(3,3)*C(3,3)/(w1*A(3,3)^2);
C4 = Q(3,3)*k2^2/(2*w1*(A(3,3))^2);
h3 = 3*(L^2)*(Fv+Kt*Kb*(n^2)/Ra)/r^2;
k = (h3^2 - 3*(L^2)*h3*Kt*Kb*(n^2)/(Ra*(r^2)))/(Jc + (3*Jw*L^2)/r^2)^2;
tau_w = 1/(k^0.5);
h2 = sinh(tf/tau_w);
h1 = 2*(1-cosh(tf/tau_w)) + h2;

% For differential equation of x_dot
R = [cos(phi(1)), -sin(phi(1)), 0; sin(phi(1)), cos(phi(1)), 0; 0, 0, 1];
R_dot = [-sin(phi(1)), -cos(phi(1)), 0; cos(phi(1)), -sin(phi(1)), 0; 0, 0, 0];
C1 = (C(1,1)/A(1,1))^2 - k2*w2*Q(1,1)*C(1,1)/(w1*A(1,1)^2);
C2 = (k2^2)*Q(1,1)/(2*w1*A(1,1)^2);
