function D = D(gamma, order, dt)
H = Bezier.H(order, dt);
z_0 = [1 zeros(1,order)];
z_T = [zeros(1,order) 1];
D = [];
for i = 1:gamma
    D = [D; z_0*H^(i-1)];
end
for i = 1:gamma
    D = [D; z_T*H^(i-1)];
end
end