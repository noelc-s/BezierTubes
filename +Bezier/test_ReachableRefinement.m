%% Description
% A script to test compare three things:
% - N step reachable over dt
% - 1 step reachable over N*dt
% - 1 step reachable with N step B-spline over N*dt

% Parameters
u_max = 5;
horizon_N = 2;
dt = 1;
A_x = [0 1; 0 -1; 1 0; -1 0];
b_x = 4*[1; 1; 1; 1];

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
[M, N, Gamma, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);
clf

%% 1 step B-spline reachable over dt, different x_bar
% Bezier Matrices for N*dt
H = Bezier.H(3, dt);
D = Bezier.D(2,3,dt);
D_nT = inv(D');
Pi = Bezier.Pi(c,2,1);
Z = Bezier.Z(3, dt);
Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;

% Forward
% Dynamic bias
x1 = X0;
x_bar = [];
for i = 1:horizon_N
    x0 = x1;
    A = Df_func(x0(1),x0(2));
    if i == 1
        x1 = x0 + (dt/horizon_N)/2*(f_func(x0(1),x0(2)));
    else
        x1 = x0 + (dt/horizon_N)*(f_func(x0(1),x0(2)));
    end
    x_bar = [x_bar x1];
end
x_barf = x_bar;

xbar = x_bar;
f_xbar = f(x_bar')';
g_xbar = g(x_bar')';

Q = Bezier.Q(horizon_N, 3);
% tic
[A_in, b_in] = Bezier.F_G(A_x, b_x, H,m, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
A = A_in;
b = b_in;


Vert_f = cddmex('extreme',struct('A',[D_vec(1:2,:); A],'B',[X0;b],'lin',1:2));
Vert_f = Bezier.Poly.conv((D_vec(3:4,:)*Vert_f.V')');
% toc

% Backward
% Dynamic bias
x1 = X0;
x_bar = [];
for i = 1:horizon_N
    x0 = x1;
    if i == 1
        x1 = x0 + (-dt/horizon_N)/2*(f_func(x0(1),x0(2)));
    else
        x1 = x0 + (-dt/horizon_N)*(f_func(x0(1),x0(2)));
    end
    x_bar = [x1 x_bar ];
end
x_barb = x_bar;

xbar = x_bar;
f_xbar = f(x_bar')';
g_xbar = g(x_bar')';

Q = Bezier.Q(horizon_N, 3);
% tic
[A_in, b_in] = Bezier.F_G(A_x, b_x, H,m, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
A = A_in;
b = b_in;

Vert_b = cddmex('extreme',struct('A',[D_vec(3:4,:); A],'B',[X0;b],'lin',1:2));
Vert_b = Bezier.Poly.conv((D_vec(1:2,:)*Vert_b.V')');

if size(Vert_b,1)>2
    ind = convhull(Vert_b);
    Vert_b = Vert_b(ind,:);
end
if size(Vert_f,1)>2
    ind = convhull(Vert_f);
    Vert_f = Vert_f(ind,:);
end

Vert1 = Vert_f;
%% N step reachable over dt, different x_bar
% Bezier Matrices for dt

H = Bezier.H(3, dt/horizon_N);
D = Bezier.D(2,3,dt/horizon_N);
D_nT = inv(D');
Pi = Bezier.Pi(c,2,1);
Z = Bezier.Z(3, dt/horizon_N);
Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;
Q = Bezier.Q(1, 3);

% subplot(3,2,4)
% hold on;

% patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
% axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-0.1 b_x(3)+0.1]);

Vert_ = X0';
tic
for k = 1:horizon_N
    Vert_tmp = [];
    for i = 1:size(Vert_,1)
        IC = Vert_(i,:)';
        
        xbar = x_barf(:,k); % to use the same xbar
        f_xbar = f(xbar');
        g_xbar = g(xbar');
        [M, N, Gamma, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);
        Pi = Bezier.Pi(c,2,1);
        
        [A_in, b_in] = Bezier.F_G(A_x, b_x, H,m, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
        A = A_in;
        b = b_in;
        
        
        x0 = IC;
        color = 'g';
        Vert = cddmex('extreme',struct('A',[D_vec(1:2,:); A],'B',[x0;b],'lin',1:2));
        Vert = Bezier.Poly.conv((D_vec(3:4,:)*Vert.V')');
        
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
    patch(Vert_(:,1),Vert_(:,2),color,'facealpha',0.1);
    
end
toc
Vert4 = Vert_;

%% Choose a trajectory for the first object, see reachable set of intermediate point

x0 = X0;
x1 = Vert1(randi(size(Vert1,1),1),:)';
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order, dt);
X = [x0 x1];
P = X(:)'/D;
Xi = [P; P*H];

Z = Bezier.Z(order, dt);
tau = linspace(0,dt/2);
C = Xi*Z(tau);

H = Bezier.H(3, dt/2);
D = Bezier.D(2,3,dt/2);
D_nT = inv(D');
Pi = Bezier.Pi(c,2,1);
Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;
Q = Bezier.Q(1, 3);

x_int = C(:,end);
xbar = x_barf(:,2); % to use the same xbar
f_xbar = f(xbar');
g_xbar = g(xbar');
[M, N, Gamma, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);
Pi = Bezier.Pi(c,2,1);

[A_in, b_in] = Bezier.F_G(A_x, b_x, H,m, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
A = A_in;
b = b_in;

x0 = x_int;
color = 'g';
Vert = cddmex('extreme',struct('A',[D_vec(1:2,:); A],'B',[x0;b],'lin',1:2));
Vert = Bezier.Poly.conv((D_vec(3:4,:)*Vert.V')');

if ~isempty(Vert)
    if size(Vert,1) > 2
        ind = convhull(Vert);
        
        Vert = Vert(ind,:);
    end
end


%% Plot results
clf
hold on;
scatter(X0(1),X0(2),50,'filled','k')
patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-1 b_x(3)+1]);
patch(Vert1(:,1),Vert1(:,2),[1 0 0],'facealpha',0.03,'edgecolor',[0.7 0.1 0.1],'linewidth',2);
patch(Vert4(:,1),Vert4(:,2),[0 0 1],'facealpha',0.03,'edgecolor',[0.1 0.1 0.7],'linewidth',2);
patch(Vert(:,1),Vert(:,2),[0 1 0],'facealpha',0.03,'edgecolor',[0.1 0.7 0.1],'linewidth',2);
legend({'','',...
    '1 step B-spline reachable over N*dt, different $$\bar{x}$$',...
    'N step reachable over dt, different $$\bar{x}$$'},...
    'interpreter','latex','fontsize',15)
plot(C(1,:),C(2,:),'linewidth',2);
tau = linspace(0,dt);
C = Xi*Z(tau);
plot(C(1,:),C(2,:),'--','linewidth',2);


