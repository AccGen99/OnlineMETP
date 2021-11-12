function val = MERV(phi_f, t, tf, tau_w, h2, h1)
    Params_Calcs;
    val = phi_f*(sinh((t-tf)/tau_w) - sinh(t/tau_w) + h2)/(tau_w*h1);
end