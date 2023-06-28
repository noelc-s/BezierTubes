function H_vec_ = H_vec(H, m, order, gamma, power_of_H)
K = Bezier.K(order+1, power_of_H+1);
H_tmp = [];
for i = 0:power_of_H
    H_tmp = [H_tmp; kron((H')^i,eye(m))];
end
H_vec_ = kron(K,eye(m))*H_tmp;
end