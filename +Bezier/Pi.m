function Pi = Pi(c,m, gamma)
C = [c(1) c(1) -c(1) -c(1);
    c(2) -c(2) c(2) -c(2)];

n = m*gamma;

A = eye(n);
B = eye(m);

Pi_1 = kron(A,ones(m,1));
Pi_2 = kron(ones(n,1), B);

Pi = [C(1,1)*Pi_1 C(2,1)*Pi_2;
    C(1,2)*Pi_1 C(2,2)*Pi_2;
    C(1,3)*Pi_1 C(2,3)*Pi_2;
    C(1,4)*Pi_1 C(2,4)*Pi_2];