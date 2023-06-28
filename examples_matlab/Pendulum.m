syms x1 x2 t
x_sym = [x1 x2];
f = [x2; -sin(x1)];
g = [0; 1];
rng('default')

% Model of system dynamics to use in controllers
f_model = f;
g_model = g;

% Symbolic gradient for MPC
Df_model = [diff(f,x1) diff(f,x2)];
Dg_model = [diff(g,x1) diff(g,x2)];

% Matlab Function-ify
f_func = matlabFunction(f,'Vars',[x_sym]);
g_func = matlabFunction(g,'Vars',[x_sym]);
Df_func = matlabFunction(Df_model,'Vars',x_sym);
Dg_func = matlabFunction(Dg_model,'Vars',x_sym);

% Define outputs
y = x1;
Dy = [diff(y,x1) diff(y,x2)];
Lfy = Dy*f_model;
Lgy = Dy*g_model;
Lf2y = [diff(Lfy,x1) diff(Lfy,x2)]*f_model;
LgLfy = [diff(Lfy,x1) diff(Lfy,x2)]*g_model;
Lf2y_func = matlabFunction(Lf2y,'Vars',x_sym);
LgLfy_func = matlabFunction(LgLfy,'Vars',x_sym);

%% Parameters
u_max = 10;
dt = 1;
A_x = [1 0; -1 0; 0 1; 0 -1];
b_x = [5; 5; 5; 5];

steps = 1;
horizon_N = 1;

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

[Q,Q_combined] = Bezier.Q(horizon_N, 3);

[M, N, Gamma, c, M_og] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);

% Bezier Matrices
order = 3;
gamma = 2;
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order,dt);
D_nT = inv(D);
Pi = Bezier.Pi(c,2,1);
Z = Bezier.Z(order, dt);

%%
pos_density = 11;
vel_density = 11;

% plot_traj = false;
% earase_plots = false;

% [H, D_nT] = Poly.getBezMatrices(order, dt);

fig = figure(1);
clf;

subplot(1,2,1)
hold on;
bR = patch(0,0,'b','facealpha',0.1);
fR = patch(0,0,'g','facealpha',0.1);
% Re = patch(0,0,'k','facealpha',0.2);
% l = patch(0,0,'k','facealpha',0.2);
% if plot_traj
%     t_fR = Poly.plotReachableTraj([0 0],'k',[0 0],[],dt,H);
%     t_bR = Poly.plotReachableTraj([0 0],'k',[0 0],[],dt,H);
%     t_Re = Poly.plotReachableTraj([0 0],'k',[0 0],[],dt,H);
% end

axis([-5 5 -5.1 5.1])
axis equal
V_x = Bezier.Poly.conv(Bezier.Poly.hyp2vert(A_x, b_x));
patch(V_x(:,1),V_x(:,2),'k','facealpha',0.05,'linewidth',2);

set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)

buffer = 1;
[X,Y] = meshgrid(linspace(-b_x(2)+buffer, b_x(1)-buffer,pos_density), linspace(-b_x(4)+buffer, b_x(3)-buffer, vel_density));

Gr = digraph;
Gr = Gr.addnode(1);
Gr.Nodes.x = [0 0];
Gr.Nodes.F = {0};
Gr.Nodes.B = {0};

ind_ = 1;
Gr = Gr.addnode(numel(X));
B_ = cell(numel(X),1);
F_ = cell(numel(X),1);
for i = 1:numel(X)
 
    x_k = [X(i); Y(i)];
    scatter(x_k(1), x_k(2));
    
    % forward x_bar
    x1 = x_k;
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
    f_xbar = f(x_bar')';
    g_xbar = g(x_bar')';
    [F,G] = Bezier.F_G(A_x, b_x, H, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
    Vert_f = cddmex('extreme',struct('A',[D(1:2,:); F],'B',[x_k;G],'lin',1:2));
    Vert_f = Vert_f.V*D(3:4,:)';
    [Hred]=cddmex('reduce_h',struct('A',[D(1:2,:); F],'B',[x_k;G],'lin',1:2));
    A_f = Hred.A;
    b_f = Hred.B;
    if size(Vert_f,1)>2 && size(uniquetol(Vert_f,1e-3,'ByRows',true),1)>1
        Vert_f = Bezier.Poly.conv(Vert_f);
        if earase_plots
            fR.XData = Vert_f(:,1);
            fR.YData = Vert_f(:,2);
        else
            patch(Vert_f(:,1),Vert_f(:,2),'g','facealpha',0.05);
        end
    else
        fR.XData = 0;
        fR.YData = 0;
    end
    
    % backward x_bar
    x1 = x_k;
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
    f_xbar = f(x_bar')';
    g_xbar = g(x_bar')';
    [F,G] = Bezier.F_G(A_x, b_x, H, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
    Vert_b = cddmex('extreme',struct('A',[D(3:4,:); F],'B',[x_k;G],'lin',1:2));
    Vert_b = Vert_b.V*D(1:2,:)';
    [Hred]=cddmex('reduce_h',struct('A',[D(3:4,:); F],'B',[x_k;G],'lin',1:2));
    A_b = Hred.A;
    b_b = Hred.B;
    if size(Vert_b,1)>2 && size(uniquetol(Vert_b,1e-3,'ByRows',true),1)>1
        Vert_b = Bezier.Poly.conv(Vert_b);
        if earase_plots
            bR.XData = Vert_b(:,1);
            bR.YData = Vert_b(:,2);
        else
            patch(Vert_b(:,1),Vert_b(:,2),'b','facealpha',0.05);
        end
    else
        bR.XData = 0;
        bR.YData = 0;
    end

        drawnow;
    Gr.Nodes.x(ind_+1,:) = x_k;
    Gr.Nodes.F{ind_+1} = [A_f b_f];
    Gr.Nodes.B{ind_+1} = [A_b b_b];
    B_{ind_} = [A_b b_b];
    F_{ind_} = [A_f b_f];
    ind_ = ind_+1;
    pause
end
Gr = rmnode(Gr,1);
drawnow;

%% Create Graph
ind = 1;
s = [];
t = [];
w = [];
for i = 1:Gr.numnodes
    F = F_{i};
    for j = 1:Gr.numnodes
        B = B_{j};
        A_in = [F(:,1:end-1); B(:,1:end-1)];
        b_in = [F(:,end); B(:,end)];
        objective=[0 0];
        IN=struct('obj',objective,'A',A_in,'B',b_in);
        OUT = cddmex('solve_lp',IN);
        if OUT.how==1
            Vert = Poly.hyp2vert(A_in, b_in);
            if size(Vert,1)>2
                Vert = Poly.conv(Vert);
                s(ind) = i;
                t(ind) = j;
                w(ind) = norm(G.Nodes.x(i,:) - mean(Vert))+norm(G.Nodes.x(j,:) - mean(Vert));
                ind = ind+1;
            end
        end
    end
end
Gr = Gr.addedge(s,t,w);
%% Plan a path
% start at the origin and end at a far eq. pt.
start_n = floor(pos_density/2)*vel_density + floor(vel_density/2)+1;
% end_n = (pos_density-1)*vel_density+ floor(vel_density/2)+1;
end_n = (pos_density-1)*floor(5*vel_density/6)+ floor(vel_density/2)+1-2;
scatter(Gr.Nodes.x(start_n,1), Gr.Nodes.x(start_n,2), 50, 'g', 'filled')
scatter(Gr.Nodes.x(end_n,1), Gr.Nodes.x(end_n,2), 50, 'b', 'filled')
path = shortestpath(Gr,start_n, end_n);
plot(Gr.Nodes.x(path,1),Gr.Nodes.x(path,2),'bo-','linewidth',2);

x_nodes = Gr.Nodes.x(path,:);
P = [];
V = [];
for i = 1:size(x_nodes,1)-1
    A_in = [Gr.Nodes.F{path(i)}(:,1:end-1); Gr.Nodes.B{path(i+1)}(:,1:end-1)];
    b_in = [Gr.Nodes.F{path(i)}(:,end); Gr.Nodes.B{path(i+1)}(:,end)];
    Vert = Poly.conv(Poly.hyp2vert(A_in, b_in));
    int_pt = mean(Vert);
    [~, pos, vel] = Poly.plotTraj('r',x_nodes(i,:),int_pt,dt,H);
    P = [P pos(1:end-1)];
    V = [V vel(1:end-1)];
    [~, pos, vel] = Poly.plotTraj('r',int_pt,x_nodes(i+1,:),dt,H);
    P = [P pos(1:end-1)];
    V = [V vel(1:end-1)];
end
subplot(1,2,2); hold on;
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)
plot(P,'linewidth',2);
plot(V,'linewidth',2);