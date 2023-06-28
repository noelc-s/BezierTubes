function D = D(gamma, order, dt)
H = Bezier.H(order, dt);
z_0 = [1; zeros(order,1)];
z_T = [zeros(order,1); 1];
D = [];
for i = 1:gamma
    D = [D H^(i-1)*z_0];
end
for i = 1:gamma
    D = [D H^(i-1)*z_T];
end
end