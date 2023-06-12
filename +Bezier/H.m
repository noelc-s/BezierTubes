function H = H(order, dt)
H = (1/dt*Bezier.R(order, dt)*Bezier.S(order,dt));
end