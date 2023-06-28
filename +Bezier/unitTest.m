%% Connect boundary conditions of minimal curves
clf
dt = 1;
gamma = 2;
order = 2*gamma-1; % minimal curve
m = 1;

x0 = [0; 0];
x1 = [1; 1];
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order, dt);
X = [x0 x1];
P = X(:)'/D;
Xi = [P; P*H];

figure(1)
scatter(Xi(1,:),Xi(2,:))
Z = Bezier.Z(order, dt);
tau = linspace(0,1);
C = Xi*Z(tau);
hold on;
plot(C(1,:),C(2,:))

%% Vectorization test
% commutation matrix K
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
assert(all(vec(Xi) == H_vec*vec(P)))

%% Refine
for r = 3:10
Xi_prime = (Xi*Bezier.R_n(3,r));
Z = Bezier.Z(r, dt);

C = Xi_prime*Z(tau);
plot(Xi_prime(1,:),Xi_prime(2,:),'r--')
scatter(Xi_prime(1,:),Xi_prime(2,:),1/r*200,'filled','r')
plot(C(1,:),C(2,:))
end

%% Split Curve
dt = 1;
gamma = 2;
order = 2*gamma-1; % minimal curve

% order_og = order;
% order = order+2;
% Xi = (Bezier.R_n(order_og,order,dt)'*Xi')';
% Z = Bezier.Z(order, dt);

segments =  10;

Q = Bezier.Q(segments, order);
for i = 1:segments
   Xi_split{i} = Xi*Q{i}'; 
end

figure(1)
scatter(Xi(1,:),Xi(2,:))
Z = Bezier.Z(order, dt);
tau = linspace(0,1);
C = Xi*Z(tau);
hold on;
plot(C(:,1),C(:,2))

for j = 1:segments
scatter(Xi_split{j}(1,:),Xi_split{j}(2,:))
Z = Bezier.Z(order, dt/2);
tau = linspace(0,dt/2);
C = Xi_split{j}*Z(tau);
hold on;
plot(C(1,:),C(2,:))
end

%% 2D curve, works now
%% Connect boundary conditions of minimal curves
clf
dt = 1;
gamma = 2;
order = 2*gamma-1; % minimal curve
m = 2;

x0 = [0; 0;... % first [pos, vel]
      0; 0];   % second [pos, vel]
x1 = [1; -1;...
      1; 1];
X = [x0 x1];
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order, dt);

x = vec(X);
Delta = Bezier.Delta(order);
Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;
p = D_vec\x;
xi = H_vec*p;

P = reshape(p,m,order+1);
Xi = reshape(xi,m*gamma,order+1);


figure(1)
hold on;
Z = Bezier.Z(order, dt);
tau = linspace(0,1);
C = Xi*Z(tau);
scatter(Xi(1,:),Xi(3,:))
scatter(Xi(2,:),Xi(4,:))
plot(C(1,:),C(3,:))
plot(C(2,:),C(4,:))

%% Vectorization test
% commutation matrix K
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
assert(all(vec(Xi) == H_vec*vec(P)))


