%% Description
% A script to test compare three things:
% - N step reachable over dt
% - 1 step reachable over N*dt
% - 1 step reachable with N step B-spline over N*dt

% Parameters
u_max = 5;
horizon_N = 5;
dt = .5;
A_x = [0 1; 0 -1; 1 0; -1 0];
b_x = 2*[1; 1; 1; 1];

% Dynamics
f = @(x) -sin(x(:,1));
% f = @(x) 0*x(:,1);
g = @(x) 1+0*x(:,1);

syms t td real
x_sym = [t td];
f_model = [td; -sin(t)];
% f_model = [td; 0];
g_model = [0; 1];
Df_model = [diff(f_model,t) diff(f_model,td)];
Dg_model = [diff(g_model,t) diff(g_model,td)];
% Matlab Function-ify
f_func = matlabFunction(f_model,'Vars',[x_sym]);
g_func = matlabFunction(g_model,'Vars',[x_sym]);
Df_func = matlabFunction(Df_model,'Vars',x_sym);
Dg_func = matlabFunction(Dg_model,'Vars',x_sym);

% Reference point
x0 = [0; 2];

xbar = x0;
f_xbar = f(xbar');
g_xbar = g(xbar');
X0 = x0;
XBAR = xbar;

Lf = 1;
Lg = 1;
e_bar = 0;
K = [-1 -1];
[M, N, Gamma, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);
clf

clf
hold on;

p_f = patch(0,0,[0 1 0],'facealpha',0.03,'edgecolor',[0.1 0.7 0.1],'linewidth',2);
p_b = patch(0,0,[0 0 1],'facealpha',0.03,'edgecolor',[0.1 0.1 0.7],'linewidth',2);
xb_f = scatter(0,0,50,[0.1 0.7 0.1],'filled');
xb_b = scatter(0,0,50,[0.1 0.1 0.7],'filled');
x0_plot = scatter(X0(1),X0(2),50,'filled','k');
patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-1 b_x(3)+1]);


%% 1 step B-spline reachable over dt, different x_bar
% Bezier Matrices for N*dt
H = Bezier.H(3, horizon_N*dt);
D = Bezier.D(2,3,horizon_N*dt);
D_nT = inv(D');
Pi = Bezier.Pi(c,2,1);
Z = Bezier.Z(3, horizon_N*dt);

for tau = 0:0.02:20*pi

X0 = [cos(tau); sin(tau)];  
    
% Forward
% Dynamic bias
x1 = X0;
x_bar = [];
for i = 1:horizon_N
    x0 = x1;
    A = Df_func(x0(1),x0(2));
    if i == 1
        x1 = expm(A*dt/2)*x0;
    else
        x1 = expm(A*dt)*x0;
    end
    x_bar = [x_bar x1];
end
x_barf = x_bar;

xbar = x_bar;
f_xbar = f(x_bar')';
g_xbar = g(x_bar')';

Q = Bezier.Q(horizon_N, 3);
% tic
[A_in, b_in] = Bezier.F_G(A_x, b_x, H, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
A = A_in;
b = b_in;


Vert_f = cddmex('extreme',struct('A',[D(1:2,:); A],'B',[X0;b],'lin',1:2));
Vert_f = Vert_f.V*D(3:4,:)';
% toc

% Backward
% Dynamic bias
x1 = X0;
x_bar = [];
for i = 1:horizon_N
    x0 = x1;
    A = Df_func(x0(1),x0(2));
    if i == 1
        x1 = expm(-A*dt/2)*x0;
    else
        x1 = expm(-A*dt)*x0;
    end
    x_bar = [x1 x_bar ];
end
x_barb = x_bar;

xbar = x_bar;
f_xbar = f(x_bar')';
g_xbar = g(x_bar')';

Q = Bezier.Q(horizon_N, 3);
% tic
[A_in, b_in] = Bezier.F_G(A_x, b_x, H, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
A = A_in;
b = b_in;

Vert_b = cddmex('extreme',struct('A',[D(3:4,:); A],'B',[X0;b],'lin',1:2));
Vert_b = Vert_b.V*D(1:2,:)';

if size(Vert_b,1)>2
    ind = convhull(Vert_b);
    Vert_b = Vert_b(ind,:);
end
if size(Vert_f,1)>2
    ind = convhull(Vert_f);
    Vert_f = Vert_f(ind,:);
end
% patch(Vert(:,1),Vert(:,2),color,'facealpha',0.1);
% Vert2 = Vert;


%% Plot results
p_f.Vertices = Vert_f;
p_f.Faces = 1:size(Vert_f,1);
p_b.Vertices = Vert_b;
p_b.Faces = 1:size(Vert_b,1);
xb_f.XData = x_barf(1,:);
xb_f.YData = x_barf(2,:);
xb_b.XData = x_barb(1,:);
xb_b.YData = x_barb(2,:);
x0_plot.XData = X0(1);
x0_plot.YData = X0(2);
drawnow
end
