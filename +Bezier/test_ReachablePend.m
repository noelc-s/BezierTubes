%% Description 
% A script to test and make sure that the reachable sets are
% underapproximations of the true reachable sets. To do this, we sample
% along the boundary of the reachable set, compute trajectories, and verify
% state and input constraint satisfaction.

% Parameters
u_max = 10;
dt = 1;
A_x = [0 1; 0 -1; 1 0; -1 0];
b_x = [1; 1; 1; 1];

% Dynamics
f = @(x) -sin(x(:,2));
g = @(x) 1+0*x(:,1);
Lf = 1;
Lg = 1;
e_bar = 0;
K = [-1 -1];
[M, N, Gamma, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);

% Reference point 
x0 = [0; 0];
xbar = [0; 0];
f_xbar = f(xbar');
g_xbar = g(xbar');

% Bezier Matrices
H = Bezier.H(3, dt);
D = Bezier.D(2,3,dt);
D_nT = inv(D');
Pi = Bezier.Pi(c,2,1);
Z = Bezier.Z(3, dt);

%% Constraint Calculation
H_0 = H^0;
H_1 = H^1;
H_2 = H^2;
A_x_ = [];
b_x_ = [];
A_u_ = [];
b_u_ = [];
for m = 1:4
    I_m = zeros(1,4);
    I_m(m) = 1;
    A_tmp = [I_m*H_0'; I_m*H_1'];
    
    % input constaints
    A_u_ = [A_u_; Pi*[eye(2)*A_tmp*D_nT;...
        I_m*H_2'*D_nT]];
    b_u_ = [b_u_; 1+Pi*[xbar; f_xbar]]; % has to match numerator of c.
    
    % and state constraints
    A_x_ = [A_x_; A_x*A_tmp*D_nT];
    b_x_ = [b_x_; b_x];
    
end
A_in = [A_u_; A_x_];
b_in = [b_u_; b_x_];

clf
IC = x0;
%% Plot Reachable sets
for r = 1:2
    if r == 1
%         Forward
        A = A_in(:,3:4);
        b = b_in - A_in(:,1:2)*IC;
        x0 = IC;
        color = 'b';
    else
%         Backward
        A = A_in(:,1:2);
        b = b_in - A_in(:,3:4)*IC;
        x1 = IC;
        color = 'g';
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
        
        xi = inv(D')*[x0; x1];
        Xi = [xi'; xi'*H];
        q_d_gamma = Z(tau)*(xi'*H^2)';
        
        subplot(2,2,r)
        scatter(Xi(1,:),Xi(2,:))
        tau = linspace(0,dt);
        X_D = Z(tau)*Xi';
        plot(X_D(:,1),X_D(:,2));
        
        U = 1./(g(X_D)).*(-f(X_D) + q_d_gamma);
        subplot(2,2,r+2)
        plot(tau,U)
    end
end
end



