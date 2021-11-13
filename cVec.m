function val = cVec(B, R, lamda_dot, alpha, omega, R_dot, lamda, C, A, w2)
    val = 0*lamda_dot;
    val = (B*inv(R)*(lamda_dot + alpha + omega*R_dot*R'*lamda - C*inv(A)*lamda))/w2;
end