%% Connect boundary conditions of minimal curves
clf
dt = 1;
gamma = 2;
order = 2*gamma-1; % minimal curve
m = 1;

x0 = [0; 0];
x1 = [.45; .7];
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order, dt);
X = [x0 x1];
p = X(:)'/D;
P = [p; p*H];

figure(1)
scatter(P(1,:),P(2,:),'filled')
Z = Bezier.Z(order, dt);
tau = linspace(0,1);
B = P*Z(tau);
hold on;
plot(B(1,:),B(2,:))

%% Vectorization test
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
assert(all(vec(P) == H_vec*vec(p)))

%% Refine
for r = 3:4
P_prime = (P*Bezier.R_n(3,r));
Z = Bezier.Z(r, dt);

B = P_prime*Z(tau);
plot(P_prime(1,:),P_prime(2,:),'r--')
scatter(P_prime(1,:),P_prime(2,:),1/r*200,'filled','r')
plot(B(1,:),B(2,:))
end

%% Split Curve
dt = 1;
gamma = 2;
order = 2*gamma-1; % minimal curve

% order_og = order;
% order = order+2;
% P = (Bezier.R_n(order_og,order,dt)'*P')';
% Z = Bezier.Z(order, dt);

segments =  5;

Q = Bezier.Q(segments, order);
for i = 1:segments
   P_split{i} = P*Q{i}; 
end
for i = 1:segments
   P_split_vec{i} = kron(Q{i}',kron(eye(m),eye(gamma)))*P(:); 
end


figure(1)
scatter(P(1,:),P(2,:))
Z = Bezier.Z(order, dt);
tau = linspace(0,1);
B = P*Z(tau);
hold on;
plot(B(:,1),B(:,2))

for j = 1:segments
scatter(P_split{j}(1,:),P_split{j}(2,:))
P_vec = reshape(P_split_vec{j},2,[]);
scatter(P_vec(1,:),P_vec(2,:),'filled')
Z = Bezier.Z(order, 1);
tau = linspace(0,1);
B = P_split{j}*Z(tau);
hold on;
plot(B(1,:),B(2,:))
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
x1 = [1; 2;...
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

p = reshape(p,m,order+1);
P = reshape(xi,m*gamma,order+1);


figure(1)
hold on;
Z = Bezier.Z(order, dt);
tau = linspace(0,1);
B = P*Z(tau);
scatter(P(1,:),P(3,:))
scatter(P(2,:),P(4,:))
plot(B(1,:),B(3,:))
plot(B(2,:),B(4,:))

%% Vectorization test
assert(all(xi == H_vec*p(:)))
assert(all(Delta_vec*xi == x))
assert(all(D_vec*p(:) == x))

%% Split Curve

segments =  10;

Q = Bezier.Q(segments, order);
for i = 1:segments
   P_split{i} = P*Q{i}; 
end

for j = 1:segments
scatter(P_split{j}(1,:),P_split{j}(3,:))
scatter(P_split{j}(2,:),P_split{j}(4,:))
Z = Bezier.Z(order, 1);
tau = linspace(0,1);
B = P_split{j}*Z(tau);
hold on;
plot(B(1,:),B(3,:))
plot(B(2,:),B(4,:))
end


%% minimize curvature of non-minimal curves
% clf
dt = 1;
gamma = 2;
order = 10;
m = 1;

x0 = [0; 0];
x1 = [1; 1];
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order, dt);
X = [x0 x1];
% Q = eye(order+1);  % same as least squares: p = X(:)'/D;
Diff = diag(-1*ones(order,1),-1)+diag([0 ones(1,order)]);
Q = Diff'*Diff;
Q = Q'*Q;
Q = H*Diff'*Diff*H';
A_eq = D';
b_eq = X(:);
p_v = [Q A_eq'; A_eq zeros(size(A_eq,1))] \ [zeros(size(Q,1),1); b_eq];
p = p_v(1:order+1)';
P = [p; p*H];

figure(1)
scatter(P(1,:),P(2,:))
Z = Bezier.Z(order, dt);
tau = linspace(0,1);
B = P*Z(tau);
hold on;
plot(B(1,:),B(2,:))