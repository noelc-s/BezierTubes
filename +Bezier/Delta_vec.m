function Delta_vec_ = Delta_vec(m, order, gamma)
K = Bezier.K(m,gamma);
Delta = Bezier.Delta(order);
Delta_vec_ = kron(Delta',K);
end