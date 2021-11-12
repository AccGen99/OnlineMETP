clc
clear
Params_Calcs;

t = linspace(0, tf, 1000);
phi_dot = 0*t;  % Array vector, phi_dot(i) corresponds to time(i)
range = 1:length(t);

for i=range
    phi_dot(i) = MERV(phi_f, t(i), tf, tau_w, h2, h1);
end
