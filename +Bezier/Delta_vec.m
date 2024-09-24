function Delta_vec_ = Delta_vec(m, order, gamma)
Delta = Bezier.Delta(order);
Delta_vec_ = kron(Delta',eye(gamma*m));
end