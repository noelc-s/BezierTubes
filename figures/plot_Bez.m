%% Description

% Parameters
u_max = 3;
dt = 1;
A_x = [1 0; -1 0; 0 1; 0 -1];
b_x = [1;1;1;1];

steps = 1;

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
% xbar = [0; 0];
% f_xbar = f(xbar');
% g_xbar = 1./g(xbar'); % This is g_inverse

syms t td real
x_sym = [t td];
f_model = [td; -gf/l*sin(t)];
g_model = [0; 1/(m*l^2)];
Df_model = [diff(f_model,t) diff(f_model,td)];
Dg_model = [diff(g_model,t) diff(g_model,td)];
% Matlab Function-ify
f_func = matlabFunction(f_model,'Vars',[x_sym]);
g_func = matlabFunction(g_model,'Vars',[x_sym]);
Df_func = matlabFunction(Df_model,'Vars',x_sym);
Dg_func = matlabFunction(Dg_model,'Vars',x_sym);

% Bezier Matrices
order = 3;
gamma = 2;
[Q,Q_combined] = Bezier.Q(steps, 3);
H = Bezier.H(3, dt);
Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;

% Dynamic bias
x1 = x0;
xbar = [];
for i = 1:steps
    x0 = x1;
    A = Df_func(x0(1),x0(2));
    if i == 1
        x1 = expm(A*(dt/steps)/2)*x0;
    else
        x1 = expm(A*dt/steps)*x0;
    end
    xbar = [xbar x1];
end
f_xbar = f(xbar')';
g_xbar = 1./g(xbar')'; % This is g_inverse

%% U_MAX
subaxis(3,3,1, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
axis off
hold on;
patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-0.1 b_x(3)+0.1]);
U_MAX = linspace(10,.1,10);

c1 = [0 1 0];
c2 = [0 0 1];
c1_f = [0.1 0.7 0.1];
c2_f = [0.1 0.1 0.7];

for j = 1:length(U_MAX)
    u_max_ = U_MAX(j);
    [A, b] = Bezier.F_G(A_x, b_x, H,m, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max_);
    % Forward
    Vert_f = cddmex('extreme',struct('A',[D_vec(1:2,:); A],'B',[x0;b],'lin',1:2));
    Vert_f = Bezier.Poly.conv((D_vec(3:4,:)*Vert_f.V')');
    if ~isempty(Vert_f)
        ind = convhull(Vert_f);
        Vert_f = Vert_f(ind,:);
    end
    patch(Vert_f(:,1),Vert_f(:,2),c1+j/length(U_MAX)*(c2-c1),'facealpha',0.01,'edgecolor',c1_f+j/length(U_MAX)*(c2_f-c1_f),'linewidth',2);
end

%% Time Horizon
subaxis(3,3,4, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
axis off
hold on;
patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-0.1 b_x(3)+0.1]);
DT = linspace(5, 0.3, 10);

c1 = [0 1 0];
c2 = [0 0 1];
c1_f = [0.1 0.7 0.1];
c2_f = [0.1 0.1 0.7];

for j = 1:length(DT)
    dt_ = DT(j);
    tic
    H = Bezier.H(3, dt_);
    Delta_vec = Bezier.Delta_vec(m, order, gamma);
    H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
    D_vec = Delta_vec*H_vec;
    
    [A, b] = Bezier.F_G(A_x, b_x, H,m, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
    % Forward
    Vert_f = cddmex('extreme',struct('A',[D_vec(1:2,:); A],'B',[x0;b],'lin',1:2));
    Vert_f = Bezier.Poly.conv((D_vec(3:4,:)*Vert_f.V')');
    if ~isempty(Vert_f)
        ind = convhull(Vert_f);
        Vert_f = Vert_f(ind,:);
    end
    patch(Vert_f(:,1),Vert_f(:,2),c1+j/length(DT)*(c2-c1),'facealpha',0.01,'edgecolor',c1_f+j/length(DT)*(c2_f-c1_f),'linewidth',2);
end

%% Initial condition
subaxis(3,3,7, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
axis off
hold on;

hold on;
patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-0.1 b_x(3)+0.1]);
U_MAX = linspace(10,.1,20);

c1 = [0 1 0];
c2 = [0 0 1];
c1_f = [0.1 0.7 0.1];
c2_f = [0.1 0.1 0.7];

num = 2;
X0 = [linspace(0,0,num); linspace(0.8,-0.8,num)];

H = Bezier.H(3, dt);
Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;


for j = 1:size(X0,2)
    
    x0 = X0(:,j);
    [A, b] = Bezier.F_G(A_x, b_x, H,m, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
    % Forward
    Vert_f = cddmex('extreme',struct('A',[D_vec(1:2,:); A],'B',[x0;b],'lin',1:2));
    Vert_f = Bezier.Poly.conv((D_vec(3:4,:)*Vert_f.V')');
    if ~isempty(Vert_f)
        ind = convhull(Vert_f);
        Vert_f = Vert_f(ind,:);
    end
    patch(Vert_f(:,1),Vert_f(:,2),c1+j/num*(c2-c1),'facealpha',0.01,'edgecolor',c1_f+j/num*(c2_f-c1_f),'linewidth',2);
    scatter(x0(1),x0(2),50,c1+j/num*(c2-c1),'filled');
end
