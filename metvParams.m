sysParams;
M = diag([m, m, J]);
B = [0, 1, L
    -0.8660, -0.5, L
    0.8660, -0.5, L];
A = M + (Jw/r^2)*(B')*B;
Q = (B')*B;
R = [cos(psi), -sin(psi), 0
    sin(psi), cos(psi), 0
    0, 0, 1];

a11 = A(1,1);
q11 = Q(1,1);
c11 = 1.5*(Fv+Kt*Kb*(n^2)/Ra)/r^2;
k2 = Vs*Kt*n/(Ra*r);
w1 = (Vs^2)/Ra;
w2 = Kb*n*Vs/(Ra*r);

C1 = (c11/a11)^2 - k2*w2*q11*c11/(w1*(a11^2));
C2 = (k2^2)*q11/(2*w1*(a11^2));