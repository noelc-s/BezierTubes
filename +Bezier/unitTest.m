%% Connect boundary conditions of minimal curves
clf
dt = 1;
gamma = 2;
order = 2*gamma-1; % minimal curve

x0 = [0 0];
x1 = [1 1];
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order, dt);
xi = inv(D')*[x0'; x1'];
Xi = [xi'; xi'*H];

figure(1)
scatter(Xi(1,:),Xi(2,:))
Z = Bezier.Z(order, dt);
tau = linspace(0,1);
C = Z(tau)*Xi';
hold on;
plot(C(:,1),C(:,2))

%% Refine
for r = 3:10
Xi_prime = (Bezier.R_n(3,r,dt)'*Xi')';
Z = Bezier.Z(r, dt);

C = Z(tau)*Xi_prime';
plot(Xi_prime(1,:),Xi_prime(2,:),'r--')
scatter(Xi_prime(1,:),Xi_prime(2,:),1/r*200,'filled','r')
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
xi = inv(D')*[x0'; x1'];
Xi = [xi'; xi'*H];

% order_og = order;
% order = order+2;
% Xi = (Bezier.R_n(order_og,order,dt)'*Xi')';
% Z = Bezier.Z(order, dt);

segments =  10;

Q = Bezier.Q(segments, order, dt);
for i = 1:segments
   Xi_split{i} = (Q{i}*Xi')'; 
end

figure(1)
scatter(Xi(1,:),Xi(2,:))
Z = Bezier.Z(order, dt);
tau = linspace(0,1);
C = Z(tau)*Xi';
hold on;
plot(C(:,1),C(:,2))

for j = 1:segments
scatter(Xi_split{j}(1,:),Xi_split{j}(2,:))
Z = Bezier.Z(order, dt/2);
tau = linspace(0,dt/2);
C = Z(tau)*Xi_split{j}';
hold on;
plot(C(:,1),C(:,2))
end

%% Higher order curve with same constrol points as Bspline
%%% Conclusion: It will not go through the same number of points, but is
%%% not so bad in the end
%%%
%%% Question: How do we metrize the difference between these curves?? They
%%% have the same control points
%%% Answer: have to split the resulting curve up and compare those
%%% segments.
clf
subplot(1,2,1)
hold on

syms t td real
x_sym = [t td];
f = [td; -sin(t)];
g = [0; 1];
Lg = 1;
rng('default')
% Model of system dynamics to use in controllers
f_model = f;
g_model = g;
% Symbolic gradient for MPC
Df_model = [diff(f,t) diff(f,td)];
Dg_model = [diff(g,t) diff(g,td)];
% Matlab Function-ify
f_func = matlabFunction(f,'Vars',[x_sym]);
g_func = matlabFunction(g,'Vars',[x_sym]);
Df_func = matlabFunction(Df_model,'Vars',x_sym);
Dg_func = matlabFunction(Dg_model,'Vars',x_sym);

dt = .1;
iter = 10; % must be divisible by refinement order at end
gamma = 2;
order = 2*gamma-1; % minimal curve

% Produce x_bar
C_true = [];
x0 = [0 1]';
A = Df_func(x0(1),x0(2));
x1 = expm(A*dt)*x0;
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order, dt);
xi = inv(D')*[x0; x1];
Xi = [xi'; xi'*H];
figure(1)
scatter(Xi(1,:),Xi(2,:),30,'b','filled')
Z = Bezier.Z(order, dt);
tau = linspace(0,dt);
C = Z(tau)*Xi';
hold on;
plot(C(:,1),C(:,2),'b','linewidth',2)
XI = Xi;
C_true = [C_true; C];
for i = 1:iter-1
x0 = x1;
A = Df_func(x0(1),x0(2));
x1 = expm(A*dt)*x0;
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order, dt);
xi = inv(D')*[x0; x1];
Xi = [xi'; xi'*H];
figure(1)
scatter(Xi(1,:),Xi(2,:),30,'filled','b')
Z = Bezier.Z(order, dt);
tau = linspace(0,dt);
C = Z(tau)*Xi';
hold on;
plot(C(:,1),C(:,2),'b','linewidth',2)
XI = [XI Xi];
C_true = [C_true; C(2:end,:)];
drawnow
end
% Z = Bezier.Z(size(XI,2)-1, dt);
% C = Z(tau)*XI';
% hold on;
% plot(C(:,1),C(:,2),'g','linewidth',2)
% plot(XI(1,:),XI(2,:),'b--','linewidth',2)
% C_orig = C;

H = Bezier.H(order, (iter)*dt);
D = Bezier.D(gamma,order, (iter)*dt);
xi = inv(D')*[XI(:,1); XI(:,end)];
Xi = [xi'; xi'*H];
scatter(Xi(1,:),Xi(2,:),30,'filled','r')
plot(Xi(1,:),Xi(2,:),'r--','linewidth',2)
Z = Bezier.Z(order, (iter)*dt);
tau = linspace(0,(iter)*dt);
C = Z(tau)*Xi';
C_fit = C;
hold on;
plot(C(:,1),C(:,2),'r','linewidth',2)
[~,Q] = Bezier.Q(iter, order);
Xi_prime = (Q*Xi')';
r = 2*iter*gamma-1;
dt_prime = iter*dt;
Xi_minimal = Xi;
Z = Bezier.Z(r, dt_prime);
C = Z(tau)*Xi_prime';
% plot(Xi_prime(1,:),Xi_prime(2,:),'r--')
scatter(Xi_prime(1,:),Xi_prime(2,:),30,'filled','r')
% plot(C(:,1),C(:,2))
% C_fit = C;

% % Fit an arbitrary number of constraints
% 
% % This is for equality constraint on control point
% splits = 2;
% X = XI(:,1:(size(XI,2)/splits-1):end);
% X = [X(:,1:end-1) XI(:,end)];
% 
% %This is for equality constraint on curve associated with control points
% % X = C_orig([1 50 100],:)';
% 
% gamma = 2;
% constraints = size(X,2);
% order = constraints*gamma-1;
% H = Bezier.H(order, (iter)*dt);
% Z = Bezier.Z(order, (iter)*dt);
% D = [];
% for k = 1:constraints
% for i = 1:gamma
%     D = [D H^(i-1)*Z((k-1)/(constraints-1)*(iter)*dt)'];
% end
% end
% xi = inv(D')*X(:);
% Xi = [xi'; xi'*H];
% figure(1)
% scatter(Xi(1,:),Xi(2,:),30,'filled','c')
% tau = linspace(0,(iter)*dt);
% C = Z(tau)*Xi';
% hold on;
% plot(C(:,1),C(:,2),'c','linewidth',2)
% plot(Xi(1,:),Xi(2,:),'c--','linewidth',2)
% r = 2*iter*gamma-1;
% dt_prime = iter*dt;
% Xi_prime_split = (Bezier.R_n(order,r,1)'*Xi')';
% clear Z
% for i = 0:r
%     Z{i+1} = @(t) nchoosek(r, i)*(t./dt_prime).^i.*(1-t./dt_prime).^(r-i);
% end
% Z = @(t) cell2mat(arrayfun(@(idx) Z{idx}(t)', 1:r+1, 'uniform',0));
% C = Z(tau)*Xi_prime_split';
% plot(Xi_prime_split(1,:),Xi_prime_split(2,:),'c--')
% scatter(Xi_prime_split(1,:),Xi_prime_split(2,:),30,'filled','c')
% plot(C(:,1),C(:,2))
% C_split = C;


% Error bounds
subplot(1,2,2)
plot(C_true(1:iter:end,:))
hold on
plot(C_fit,'r')
plot(vecnorm(C_fit'-C_true(1:iter:end,:)'),'c--')
Bez_bound = max(vecnorm(XI-Xi_prime));
line([0 100], [Bez_bound Bez_bound],'color','c','linewidth',2);
Ctrl_bound = Lg*max(vecnorm(Xi_prime(:,order+1:order+1:end) - Xi_prime(:,1:order+1:end-(order)))) + Lg*max(vecnorm(XI-Xi_prime));
line([0 100], [Ctrl_bound Ctrl_bound],'color','r','linewidth',2);

% plot(C_split,'c')
% plot(vecnorm(C_split'-C_orig'),'c--')
% Ctrl_bound = max(vecnorm(XI-Xi_prime_split));
% line([0 100], [Ctrl_bound Ctrl_bound],'color','c','linewidth',2);

%%%prior bound
prior_bound = Lg*max(vecnorm(Xi_minimal - XI(:,1)));
line([0 100], [prior_bound prior_bound],'color','k','linewidth',2);

