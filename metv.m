function dydt = metv(t, y, alpha_x)
    metvParams;
    w = merv(t);
    dydt(1) = y(2);
    dydt(2) = ((w^2)+C1)*y(1)+C2*alpha_x;
end