function val = maxAcc(k2, R, A, B, u_max)
    % Max acceleration vector for t = 0
    x_dd = k2*inv(R')*inv(A)*B'*u_max;
    val = x_dd(1);
end