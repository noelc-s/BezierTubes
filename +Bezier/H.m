function H = H(order, dt)
H = (1/dt*Bezier.S(order,dt)'*Bezier.R(order, dt)');
end