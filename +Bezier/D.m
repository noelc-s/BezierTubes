function D = D(gamma, order, dt)
H = Bezier.H(order, dt);
z_0 = [1 zeros(1,order)]';
z_T = [zeros(1,order) 1]';
D = [];
for i = 1:gamma
    D = [D H^(i-1)*z_0];
end
for i = 1:gamma
    D = [D H^(i-1)*z_T];
end
end