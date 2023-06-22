%% Connect boundary conditions of minimal curves
clf
dt = 1;
gamma = 2;
order = 2*gamma-1; % minimal curve

x0 = [0 0];
x1 = [1 1];
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order, dt);
xi = D\[x0'; x1'];
Xi = [xi H*xi];

figure(1)
scatter(Xi(:,1),Xi(:,2))
Z = Bezier.Z(order, dt);
tau = linspace(0,1);
C = Z(tau)'*Xi;
hold on;
plot(C(:,1),C(:,2))

%% Refine
for r = 3:10
Xi_prime = (Bezier.R_n(3,r)*Xi);
Z = Bezier.Z(r, dt);

C = Z(tau)'*Xi_prime;
plot(Xi_prime(:,1),Xi_prime(:,2),'r--')
scatter(Xi_prime(:,1),Xi_prime(:,2),1/r*200,'filled','r')
plot(C(:,1),C(:,2))
end

%% Split Curve
dt = 1;
gamma = 2;
order = 2*gamma-1; % minimal curve

x0 = [0 0];
x1 = [1 1];
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order, dt);
Z = Bezier.Z(order, dt);
xi = inv(D)*[x0'; x1'];
Xi = [xi H*xi];

% order_og = order;
% order = order+2;
% Xi = (Bezier.R_n(order_og,order,dt)'*Xi')';
% Z = Bezier.Z(order, dt);

segments =  10;

Q = Bezier.Q(segments, order);
for i = 1:segments
   Xi_split{i} = Q{i}*Xi; 
end

figure(1)
scatter(Xi(1,:),Xi(2,:))
Z = Bezier.Z(order, dt);
tau = linspace(0,1);
C = Z(tau)'*Xi;
hold on;
plot(C(:,1),C(:,2))

for j = 1:segments
scatter(Xi_split{j}(:,1),Xi_split{j}(:,2))
Z = Bezier.Z(order, dt/2);
tau = linspace(0,dt/2);
C = Z(tau)'*Xi_split{j};
hold on;
plot(C(:,1),C(:,2))
end

%% 2D curve
%% Connect boundary conditions of minimal curves
clf
dt = 1;
gamma = 2;
order = 2*gamma-1; % minimal curve

x0 = [0 0;...
      0 0];
x1 = [1 1;...
      1 1];
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order, dt);
xi = D\[x0(:,1) x0(:,2) x1(:,1) x1(:,2)]';


Xi = [xi H*xi];

figure(1)
scatter(Xi(:,1),Xi(:,2))
Z = Bezier.Z(order, dt);
tau = linspace(0,1);
C = Z(tau)'*Xi;
hold on;
plot(C(:,1),C(:,2))



