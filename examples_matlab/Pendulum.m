clear;
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
u_max = 9;
dt = 1;
A_x = [1 0; -1 0; 0 1; 0 -1];
b_x = 5*[1;1;1;1];

% steps = 1;
horizon_N = 10;

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
Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;

%%
% pos_density = 9;
% vel_density = 9;
pos_bounds = [-5 5];
vel_bounds = [-5 5];
% pos_bounds = [pi pi];
% vel_bounds = [0 0];
delta_pos = .5;
delta_vel = .5;

plot_traj = false;
earase_plots = true;
plot_anything = false;

% [H, D_nT] = Poly.getBezMatrices(order, dt);

fig = figure(1);
clf;

subplot(2,2,1)
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
[X,Y] = meshgrid([pos_bounds(1):delta_pos:pos_bounds(2) pi],...
    [vel_bounds(1):delta_vel:vel_bounds(2) 0]);

Gr = digraph;
Gr = Gr.addnode(1);
Gr.Nodes.x = [0 0];
Gr.Nodes.F = {0};
Gr.Nodes.B = {0};

ind_ = 1;
Gr = Gr.addnode(numel(X));
B_ = cell(numel(X),1);
F_ = cell(numel(X),1);
X_K = cell(numel(X),1);
for i = 1:numel(X)
    disp(100*i/numel(X));
 
    x_k = [X(i); Y(i)];
    if plot_anything
    scatter(x_k(1), x_k(2));
    end
    
    % forward x_bar
    x1 = x_k;
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
    f_xbar = f(x_bar')';
    g_xbar = g(x_bar')';
    [F,G] = Bezier.F_G(A_x, b_x, H,m, x_bar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
    Vert_f = cddmex('extreme',struct('A',[D_vec(1:2,:); F],'B',[x_k;G],'lin',1:2));
    Vert_f = (D_vec(3:4,:)*Vert_f.V')';
    if size(Vert_f,1)>2 && size(uniquetol(Vert_f,1e-3,'ByRows',true),1)>1
        Vert_f = Bezier.Poly.conv(Vert_f);
        V=struct('V',[Vert_f]);
        Hyp=cddmex('hull',V);
        A_f = Hyp.A;
        b_f = Hyp.B; 
        if plot_anything && earase_plots
            fR.XData = Vert_f(:,1);
            fR.YData = Vert_f(:,2);
        elseif ~earase_plots
            patch(Vert_f(:,1),Vert_f(:,2),'g','facealpha',0.05);
        end
    else
        fR.XData = 0;
        fR.YData = 0;
        A_f = [];
        b_f = [];
    end
    
    % backward x_bar
    x1 = x_k;
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
    f_xbar = f(x_bar')';
    g_xbar = g(x_bar')';
    [F,G] = Bezier.F_G(A_x, b_x, H,m, x_bar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
    Vert_b = cddmex('extreme',struct('A',[D_vec(3:4,:); F],'B',[x_k;G],'lin',1:2));
    Vert_b = (D_vec(1:2,:)*Vert_b.V')';
    if size(Vert_b,1)>2 && size(uniquetol(Vert_b,1e-3,'ByRows',true),1)>1
        Vert_b = Bezier.Poly.conv(Vert_b);
        V=struct('V',[Vert_b]);
        Hyp=cddmex('hull',V);
        A_b = Hyp.A;
        b_b = Hyp.B; 
        if norm(A_b) > 1e5
            disp("uh oh");
        end
        if plot_anything && earase_plots
            bR.XData = Vert_b(:,1);
            bR.YData = Vert_b(:,2);
        elseif ~earase_plots
            patch(Vert_b(:,1),Vert_b(:,2),'b','facealpha',0.05);
        end
    else
        bR.XData = 0;
        bR.YData = 0;
        A_b = [];
        b_b = [];
    end

        drawnow;
    Gr.Nodes.x(ind_+1,:) = x_k;
    Gr.Nodes.F{ind_+1} = [A_f b_f];
    Gr.Nodes.B{ind_+1} = [A_b b_b];
    X_K{ind_} = x_k;
    B_{ind_} = [A_b b_b];
    F_{ind_} = [A_f b_f];
    ind_ = ind_+1;
%     pause
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
    if ~isempty(F)
        for j = 1:Gr.numnodes
            
            B = B_{j};
            if ~isempty(B)
                A_in = [F(:,1:end-1); B(:,1:end-1)];
                b_in = [F(:,end); B(:,end)];
                objective=[0 0];
                %         cddmex('extreme',struct('A',[D_vec(1:2,:); F],'B',[x_k;G],);
%                 IN= struct('obj',objective,'A',[D_vec;A_in],'B',[X_K{i};X_K{j};b_in],'lin',1:4);
                IN= struct('obj',objective,'A',[A_in],'B',[b_in]);
                OUT = cddmex('solve_lp',IN);
                if OUT.how==1
                    Vert = Bezier.Poly.hyp2vert(A_in, b_in);
                    if size(Vert,1)>2
                        Vert = Bezier.Poly.conv(Vert);
                        s(ind) = i;
                        t(ind) = j;
                        if plot_anything && plot_traj
                            line([Gr.Nodes.x(i,1) Gr.Nodes.x(j,1)],[Gr.Nodes.x(i,2) Gr.Nodes.x(j,2)])
                        end
                        w(ind) = norm(Gr.Nodes.x(i,:) - mean(Vert))+norm(Gr.Nodes.x(j,:) - mean(Vert));
%                         w(ind) = 1;
                        ind = ind+1;
                    end
                end
            end
        end
    end
end
Gr = Gr.addedge(s,t,w);
%% Plan a path
x0_des = [0 0];
x1_des = [pi 0];
[~,start_n] = min(vecnorm((Gr.Nodes.x-x0_des)',2));
[~,end_n] = min(vecnorm((Gr.Nodes.x-x1_des)',2));
% start at the origin and end at a far eq. pt.
% start_n = floor(pos_density/2)*vel_density + floor(vel_density/2)+1;
% end_n = (pos_density-1)*vel_density+ floor(vel_density/2)+1;
% end_n = (pos_density-1)*floor(5/6*vel_density)+ floor(vel_density/2)+1;
subplot(2,2,1)
scatter(Gr.Nodes.x(start_n,1), Gr.Nodes.x(start_n,2), 50, 'g', 'filled')
scatter(Gr.Nodes.x(end_n,1), Gr.Nodes.x(end_n,2), 50, 'b', 'filled')
path = shortestpath(Gr,start_n, end_n);
plot(Gr.Nodes.x(path,1),Gr.Nodes.x(path,2),'bo-','linewidth',2);

x_nodes = Gr.Nodes.x(path,:);
Pos = [];
V = [];
T = [0];
for i = 1:size(x_nodes,1)-1
    A_in = [Gr.Nodes.F{path(i)}(:,1:end-1); Gr.Nodes.B{path(i+1)}(:,1:end-1)];
    b_in = [Gr.Nodes.F{path(i)}(:,end); Gr.Nodes.B{path(i+1)}(:,end)];
    Vert = Bezier.Poly.conv(Bezier.Poly.hyp2vert(A_in, b_in));
    int_pt = mean(Vert);
    
    x0 = x_nodes(i,:)';
    x1 = int_pt';
    X = [x0 x1];
    P = X(:)'/D;
    Xi = [P; P*H];
    tau = linspace(0,dt);
    C = Xi*Z(tau);
    Pos = [Pos C(1,:)];
    V = [V C(2,:)];
    T = [T tau+T(end)];
    
    x0 = int_pt';
    x1 =  x_nodes(i+1,:)';
    X = [x0 x1];
    P = X(:)'/D;
    Xi = [P; P*H];
    tau = linspace(0,dt);
    C = Xi*Z(tau);
    Pos = [Pos C(1,:)];
    V = [V C(2,:)];
    T = [T tau+T(end)];
end
T = T(2:end);
subplot(2,2,[3 4]); hold on;
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)
plot(T,Pos,'linewidth',2);
plot(T,V,'linewidth',2);

%% Animate
subplot(2,2,2)
set(gcf,'renderer','painters')

% axis off
hold on;
P1_out = line([0 l],[0 l],'color','k','marker','o','markerSize',10,'linewidth',20);
P1_in = line([0 l],[0 l],'color',[0 0.4470 0.7410],'marker','o','markerSize',10,'linewidth',15);
axis([-1.5 1.5 -1.5 1.5])
axis equal
buff = 5;
set(gca,'FontSize',17)
set(gca,'linewidth',2)

tic
cont = true;
ind = 1;
while cont
    set(P1_out,'XData',[0 l*sin(Pos(ind))])
    set(P1_out,'YData',[0 -l*cos(Pos(ind))])
    set(P1_in,'XData',[0 l*sin(Pos(ind))])
    set(P1_in,'YData',[0 -l*cos(Pos(ind))])
    drawnow
    time = toc;
    axis([-1.5 1.5 -1.5 1.5])
        
    ind = find(T>time,1,'first');
    if isempty(ind)
        cont = false;
    end
end