%% Description
% A script to test compare three things:
% - N step reachable over dt
% - 1 step reachable over N*dt
% - 1 step reachable with N step B-spline over N*dt
clear;clc;

% Parameters
u_max =1;
horizon_N = 4;
dt = .25;
A_x = [1 0; -1 0; 0 1; 0 -1];
b_x = [0.3; 0.3; 0.3; 0.3];

% Dynamics
f = @(x) -sin(x(:,2));
g = @(x) 1+0*x(:,1);

syms t td real
x_sym = [t td];
f_model = [td; -sin(t)];
g_model = [0; 1];
Df_model = [diff(f_model,t) diff(f_model,td)];
Dg_model = [diff(g_model,t) diff(g_model,td)];
% Matlab Function-ify
f_func = matlabFunction(f_model,'Vars',[x_sym]);
g_func = matlabFunction(g_model,'Vars',[x_sym]);
Df_func = matlabFunction(Df_model,'Vars',x_sym);
Dg_func = matlabFunction(Dg_model,'Vars',x_sym);

% Reference point
x0 = [0; 0];

xbar = x0;
f_xbar = f(xbar');
g_xbar = g(xbar');
X0 = x0;
XBAR = xbar;

Lf = 1;
Lg = 1;
e_bar = 0;
K = [-1 -1];
clf

m=1;
gamma = 2;
order = 3;


%% 1 step B-spline reachable over N*dt
% Bezier Matrices for N*dt
H = Bezier.H(3, horizon_N*dt);
Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;

[Q,Q_combined] = Bezier.Q(1, 3);
[A_in, b_in] = Bezier.F_G(A_x, b_x, H, m, xbar, f_xbar, g_xbar, 2,Q,Lg, Lf, e_bar, K, u_max);
A = A_in;
b = b_in;

% Backward
Vert = cddmex('extreme',struct('A',[D_vec(3:4,:); A],'B',[X0;b],'lin',1:2));
Vert = Bezier.Poly.conv((D_vec(1:2,:)*Vert.V')');

if ~isempty(Vert)
    ind = convhull(Vert);
    Vert = Vert(ind,:);
end
Vert1 = Vert;

%% 1 step B-spline reachable over dt, different x_bar
H = Bezier.H(3, horizon_N*dt);
Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;

% Dynamic bias
x1 = X0;
x_bar = [];
for i = 1:horizon_N
    x0 = x1;
    if i == 1
        x1 = x0 + (dt/horizon_N)/2*(f_func(x0(1),x0(2)));
    else
        x1 = x0 + (dt/horizon_N)*(f_func(x0(1),x0(2)));
    end
    x_bar = [x_bar x1];
end

% Kinematic Bias
% x1 = X0;
% x_bar = [];
% for i = 1:horizon_N
%     x0 = x1;
%     x1 = x0 - [0; 0.1];
%     x_bar = [x_bar x1];
% end

xbar = x_bar;
f_xbar = f(x_bar')';
g_xbar = g(x_bar')';

Q = Bezier.Q(horizon_N, 3);
[A_in, b_in] = Bezier.F_G(A_x, b_x, H,m, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
A = A_in;
b = b_in;

% Backward
Vert = cddmex('extreme',struct('A',[D_vec(3:4,:); A],'B',[X0;b],'lin',1:2));
Vert = Bezier.Poly.conv((D_vec(1:2,:)*Vert.V')');

if ~isempty(Vert)
    ind = convhull(Vert);
    Vert = Vert(ind,:);
end
Vert2 = Vert;

%% N step reachable over dt, same x_bar

H = Bezier.H(3, dt);
D = Bezier.D(2,3,dt);

Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;


Vert_ = X0';

xbar = X0;
f_xbar = f(xbar');
g_xbar = g(xbar');
[M, N, Gamma, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);

for k = 1:horizon_N
    Vert_tmp = [];
    for i = 1:size(Vert_,1)
        IC = Vert_(i,:)';
        
        [F, G] = Bezier.F_G(A_x, b_x, H, m, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);

        x0 = IC;
        color = 'g';
        Vert = cddmex('extreme',struct('A',[D_vec(3:4,:); F],'B',[IC;G],'lin',1:2));
        Vert = (D_vec(1:2,:)*Vert.V')';
        
        if ~isempty(Vert)
            if size(Vert,1) > 2
                ind = convhull(Vert);
                
                Vert = Vert(ind,:);
            end
            Vert_tmp = [Vert_tmp; Vert(1:end-1,:)];
        end
    end
    Vert_ = Vert_tmp;
        if isempty(Vert_)
            Vert_ = ones(0,2);
        break
    end
    ind = convhull(Vert_);
    Vert_ = Vert_(ind,:);
    
end
Vert3 = Vert_;

%% N step reachable over dt, different x_bar

H = Bezier.H(3, dt);
Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;
Vert_ = X0';
for k = 1:horizon_N
    Vert_tmp = [];
    for i = 1:size(Vert_,1)
        IC = Vert_(i,:)';
        
        xbar = IC;
        f_xbar = f(xbar');
        g_xbar = g(xbar');
        [M, N, Gamma, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);
        
        [F, G] = Bezier.F_G(A_x, b_x, H, m, xbar, f_xbar, g_xbar, gamma,{eye(4)},Lg,Lf,e_bar,K,u_max);
        A = F;
        b = G;
        
        
        x0 = IC;
        color = 'g';
%         Vert = cddmex('extreme',struct('A',[D(1:2,:); A],'B',[x0;b],'lin',1:2));
%         Vert = uniquetol(Vert.V*D(3:4,:)',1e-3,'ByRows',1);
        
        Vert = cddmex('extreme',struct('A',[D_vec(3:4,:); F],'B',[IC;G],'lin',1:2));
        Vert = (D_vec(1:2,:)*Vert.V')';
        
        if ~isempty(Vert)
            if size(Vert,1) > 2
                ind = convhull(Vert);
                Vert = Vert(ind,:);
            end
            Vert_tmp = [Vert_tmp; Vert(1:end-1,:)];
        end
    end
    Vert_ = Vert_tmp;
    if isempty(Vert_)
        break
    end
    ind = convhull(Vert_);
    Vert_ = Vert_(ind,:);
    
    
end

[F, G] = Bezier.F_G(A_x, b_x, H, m, xbar, f_xbar, g_xbar, gamma,{eye(4)},Lg,Lf,e_bar,K,u_max);

[Hred]=cddmex('hull',struct('V',Vert_));
A = Hred.A;
b = Hred.B;

Vert = cddmex('extreme',struct('A',[[A*D_vec(1:2,1:2) zeros(size(A,1),2)]; F(end-16:end,:)],'B',[b;G(end-16:end)]));
Vert_ = (D_vec(1:2,:)*Vert.V')';
ind = convhull(Vert_);
Vert_ = Vert_(ind,:);
Vert4 = Vert_;

%% Plot results
clf
hold on;
scatter(X0(1),X0(2),50,'filled','k')
patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-.1 b_x(3)+.1]);
patch(Vert1(:,1),Vert1(:,2),[1 0 0],'facealpha',0.03,'edgecolor',[0.7 0.1 0.1],'linewidth',2);
patch(Vert2(:,1),Vert2(:,2),[0 1 0],'facealpha',0.03,'edgecolor',[0.1 0.7 0.1],'linewidth',2);
patch(Vert3(:,1),Vert3(:,2),[0 0 0],'facealpha',0.03,'edgecolor',[0 0 0],'linewidth',2);
patch(Vert4(:,1),Vert4(:,2),[0 0 1],'facealpha',0.03,'edgecolor',[0.1 0.1 0.7],'linewidth',2);
legend({'','','1 step reachable over N*dt', ...
    '1 step B-spline reachable over dt, different x_bar',...
    'N step reachable over dt, same x_bar', ...
    'N step reachable over dt, different x_bar'}, ...
    'interpreter','latex','fontsize',15)
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)
title('$u_{max} = 10$','interpreter','latex');
axis square
axis equal

%% True reachability
% f = @(x) [x(1,:); -sin(x(2,:))];
% g = @(x) 1+0*x(1,:);
% 
% u_discretization = 50;
% X_FL = [];
% X_U = [];
% tic
% x0 = X0;
% 
% tic
% for i = 1:u_discretization
%     for j = 1:u_discretization
%         %         for k = 1:u_discretization
%         %             u = @(t) t.^2*(-u_max + (i-1)/(u_discretization-1)*2*u_max)...
%         %                 +2*t.*(1-t)*(-u_max + (j-1)/(u_discretization-1)*2*u_max)...
%         %                 +(1-t).^2*(-u_max + (k-1)/(u_discretization-1)*2*u_max);
%         u = @(t) t*(-u_max + (i-1)/(u_discretization-1)*2*u_max)...
%             +(1-t)*(-u_max + (j-1)/(u_discretization-1)*2*u_max);
%         % Input
%         [t,x] = ode45(@(t,x) f(x)+g(x)*u(t),[0 horizon_N*dt], x0);
%         state_constraint_satisfaction = all(sum(A_x*x'-b_x <= 0,1) == 4);
%         input_constraint_satisfaction = 1; % by construction
%         if state_constraint_satisfaction && input_constraint_satisfaction
%             X_U = [X_U; x(end,:)];
%         end
%     end
% end
% toc
% tic
% for i = 1:u_discretization
%     for j = 1:u_discretization
%         u = @(t) t*(-u_max + (i-1)/(u_discretization-1)*2*u_max)...
%             +(1-t)*(-u_max + (j-1)/(u_discretization-1)*2*u_max);
%         % Feedback linearized
%         [t,x] = ode45(@(t,x) [0 1; 0 0]*x + [0; 1]*u(t),[0 horizon_N*dt], x0);
%         state_constraint_satisfaction = all(sum(A_x*x'-b_x <= 0,1) == 4);
%         input_constraint_satisfaction = all(abs(1./g(x').*(-[0 1]*f(x') + u(t)'))<u_max);
%         if state_constraint_satisfaction && input_constraint_satisfaction
%             X_FL = [X_FL; x(end,:)];
%         end
%         %         end
%     end
% end
% toc
% %%
% V_FL = Bezier.Poly.conv(X_FL);
% V_U= Bezier.Poly.conv(X_U);
% patch(V_FL(:,1),V_FL(:,2),[1 1 0],'facealpha',0.03,'edgecolor',[0.7 0.7 0.1],'linewidth',2);
% patch(V_U(:,1),V_U(:,2),[0 1 1],'facealpha',0.03,'edgecolor',[0.1 0.7 0.7],'linewidth',2);
% l = legend;
% l.String{end-1} = 'True Reachable Set of Linear Feedback linearized Inputs';
% l.String{end} = 'True Reachable Set of Linear Inputs';