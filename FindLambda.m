function val = FindLambda(w2, A, x_dot, k2, w1, Q, x_dot_dot, C, omega, R, R_dot)
    val = w2*(A*x_dot)/k2 - 2*w1*(A^2)*(inv(Q)*x_dot_dot) - 2*w1*A*C*(inv(Q)*x_dot)/(k2^2) + 2*w1*omega*(A^2)*inv(Q)*R*R_dot'*x_dot;
end 
