function w = merv(t)
    sysParams;
    h3 = 3*(L^2)*(Fv+Kt*Kb*(n^2)/Ra)/(r^2);
    k = ((h3^2)-3*(L^2)*h3*Kt*Kb*(n^2)/(Ra*(r^2)))/(Jc+3*Jw*((L/r)^2));
    tau_w = k^0.5;
    
    h2 = (exp(tf/tau_w) - exp(-tf/tau_w))/2;
    h1 = 2*(1-(exp(tf/tau_w) + exp(-tf/tau_w))/2) + h2;
    
    val1 = (t-tf)/tau_w;
    val2 = t/tau_w;
    w = (psi_f*(exp(val1) - exp(-val1))/2 - (exp(val2) - exp(-val2))/2 + h2)/(tau_w*h1);
%    sinh(x) = exp(x) - exp(-x) / 2
end