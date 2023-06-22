%% Description 
% A script to test and make sure that the reachable sets are
% underapproximations of the true reachable sets. To do this, we sample
% along the boundary of the reachable set, compute trajectories, and verify
% state and input constraint satisfaction.

% Parameters
u_max = 10;
dt = 1;
A_x = [1 0; -1 0; 0 1; 0 -1];
b_x = [0.3; 0.3; 0.6; 0.6];

steps = 20;

m = 1;
l = 1;
gf = 1;

% Dynamics
f = @(x) -gf/l*sin(x(:,1));
g = @(x) 1/(m*l^2)+0*x(:,1);
Lf = 1;
Lg = 1; % this is LG_inverse
e_bar = 0;
K = [-1 -1];
% Reference point 
x0 = [0; 0];
xbar = [0; 0];
f_xbar = f(xbar');
g_xbar = 1./g(xbar'); % This is g_inverse

[Q,Q_combined] = Bezier.Q(steps, 3);

[M, N, Gamma, c, M_og] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);

% Bezier Matrices
order = 3;
gamma = 2;
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order,dt);
D_nT = inv(D);
Pi = Bezier.Pi(c,2,1);
Z = Bezier.Z(order, dt);

xbar = repmat(xbar,1,steps);
f_xbar = repmat(f_xbar,1,steps);
g_xbar = repmat(g_xbar,1,steps);

%% Constraint Calculation

[F, G] = Bezier.F_G(A_x, b_x, Pi, H, xbar, f_xbar, g_xbar, gamma,Q,Lg,Lf,e_bar,K,u_max);

A_in = F*D_nT;
b_in = G;

clf
IC = x0;
%% Plot Reachable sets
for r = 1:2
    if r == 1
%         Forward
        A = A_in(:,3:4);
        b = b_in - A_in(:,1:2)*IC;
        x0 = IC;
        color = 'g';
    else
%         Backward
        A = A_in(:,1:2);
        b = b_in - A_in(:,3:4)*IC;
        x1 = IC;
        color = 'b';
    end

subplot(2,2,r)
hold on;
patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-0.1 b_x(3)+0.1]);
Vert = Bezier.Poly.hyp2vert(A,b);
ind = convhull(Vert);
Vert = Vert(ind,:);
patch(Vert(:,1),Vert(:,2),color,'facealpha',0.1);
subplot(2,2,r+2)
hold on;
line([0 dt],[1 1]*u_max)
line([0 dt],-[1 1]*u_max)
axis([0 dt -u_max-1 u_max+1]);

for i = 1:size(Vert,1)-1
    for lambda = 0:0.1:1
        if r == 1
            x1 = (lambda*Vert(i,:) + (1-lambda)*Vert(i+1,:))';
        else
            x0 = (lambda*Vert(i,:) + (1-lambda)*Vert(i+1,:))';
        end
        
        xi = inv(D)*[x0; x1];
        Xi = [xi H*xi];
        q_d_gamma = Z(tau)'*(H^2*xi);
        
        subplot(2,2,r)
        tau = linspace(0,dt);
        X_D = Z(tau)'*Xi;
        plot(X_D(:,1),X_D(:,2));
        Xi = Q_combined*Xi;
        scatter(Xi(:,1),Xi(:,2))
        
        U = 1./(g(X_D)).*(-f(X_D) + q_d_gamma);
        subplot(2,2,r+2)
        plot(tau,U)
    end
end
end



