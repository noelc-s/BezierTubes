syms x1 x2 t
x_sym = [x1 x2];
f = [x2; sin(x1)];
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

% FL
syms v
u = LgLfy\(-Lf2y + v);
FL_u = matlabFunction(u,'Vars',[x1, x2, v]);

% CLF outputs
eta = matlabFunction([y; Lfy],'Vars',x_sym);
clear x1 x2

Lf = .1;
Lg = .1;

g_xbar = 1;

M_og = 1/2*[2*Lg*Lf Lg; Lg 0];
N = [Lf*norm(g_xbar,2); norm(g_xbar,2)];
M = projectOntoSemidefiniteCone(M_og);
%%


plot_traj = false;
earase_plots = false;

pos_density = 11;
vel_density = 11;
u_max = 2;
dt = 1;

p1 = [(-N(1) + sqrt(N(1)^2+4*M(1,1)*u_max))/(2*M(1,1)) 0];
p2 = [0 (-N(2) + sqrt(N(2)^2+4*M(2,2)*u_max))/(2*M(2,2))];
c_ = [1/p1(1); 1/p2(2)];

A_x = [1 0; -1 0; 0 1; 0 -1];
b_x = [5; 5; 5; 5];

order = 3;
[H, D_nT] = Poly.getBezMatrices(order, dt);

f = figure(1);
clf;

subplot(1,2,1)
hold on;
bR = patch(0,0,'b','facealpha',0.1);
fR = patch(0,0,'g','facealpha',0.1);
Re = patch(0,0,'k','facealpha',0.2);
l = patch(0,0,'k','facealpha',0.2);
if plot_traj
    t_fR = Poly.plotReachableTraj([0 0],'k',[0 0],[],dt,H);
    t_bR = Poly.plotReachableTraj([0 0],'k',[0 0],[],dt,H);
    t_Re = Poly.plotReachableTraj([0 0],'k',[0 0],[],dt,H);
end

axis([-5 5 -5.1 5.1])
axis equal
V_x = lcon2vert(A_x,b_x);
ind_ux = convhull(V_x);
V_x = V_x(ind_ux,:);
patch(V_x(:,1),V_x(:,2),'k','facealpha',0.05,'linewidth',2);

set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)

[X,Y] = meshgrid(linspace(-b_x(2)+0.05, b_x(1)-0.05,pos_density), linspace(-b_x(4)+0.05, b_x(3)-0.05, vel_density));

G = digraph;
G = G.addnode(1);
G.Nodes.x = [0 0];
G.Nodes.F = {0};
G.Nodes.B = {0};

ind_ = 1;
G = G.addnode(numel(X));
B_ = cell(numel(X),1);
F_ = cell(numel(X),1);
for i = 1:numel(X)
 
    x_center = [X(i) Y(i)];
%     scatter(x_center(1), x_center(2));
    x_bar = x_center';
    f_xbar = f_func(x_bar(1), x_bar(2));
    g_xbar = g_func(x_bar(1), x_bar(2));

    %%% Forward Reachable
    x1 = x_center;
    A_x_ = [];
    b_x_ = [];
    A_u = [];
    b_u = [];
    for m = 1:4
        I_m = zeros(1,4);
        I_m(m) = 1;
        Ctrl_m = [I_m*H^0'; I_m*H^1']*D_nT;

        A_in = [c_(1)*Ctrl_m; c_(2)*I_m*H^2'*D_nT];
        b_in = [c_(1)*x_bar; -c_(2)*[0 1]*f_xbar];

        % forward
%         A_x_ = [A_x_; A_x*Ctrl_m(:,3:4)];
%         b_x_ = [b_x_; b_x];
%         A_u = [A_u; A_in(:,3:4); -A_in(:,3:4)];
%         b_u = [b_u; 1+b_in-A_in(:,1:2)*x0'; 1-b_in+A_in(:,1:2)*x0'];

        % backward
        A_x_ = [A_x_; A_x*Ctrl_m(:,1:2)];
        b_x_ = [b_x_; b_x];
        A_u = [A_u; A_in(:,1:2); -A_in(:,1:2)];
        b_u = [b_u; 1+b_in-A_in(:,3:4)*x0'; 1-b_in+A_in(:,3:4)*x0'];
    end
    %%% Then add constraints such that x1, x2, and x3 satisfy the constraints
    %%% too
    A_b = [A_u; A_x_];
    b_b = [b_u; b_x_];
    V = Poly.hyp2vert(A_b, b_b);
    if size(V,1)>2
        V = Poly.conv(V);
        if earase_plots
            bR.XData = V(:,1);
            bR.YData = V(:,2);
        else
            patch(V(:,1),V(:,2),'b','facealpha',0.05);
        end
    else
        bR.XData = 0;
        bR.YData = 0;
    end

    %%% Reverse Reachable
    x0 =x_center;
    A_x_ = [];
    b_x_ = [];
    A_u = [];
    b_u = [];
    for m = 1:4

        I_m = zeros(1,4);
        I_m(m) = 1;
        Ctrl_m = [I_m*H_0'; I_m*H_1']*D_nT;

        A_in = [c_(1)*Ctrl_m; c_(2)*I_m*H_2'*D_nT];
        b_in = [c_(1)*x_bar; -c_(2)*[0 1]*f_xbar];

        % forward
        A_x_ = [A_x_; A_x*Ctrl_m(:,3:4)];
        b_x_ = [b_x_; b_x];
        A_u = [A_u; A_in(:,3:4); -A_in(:,3:4)];
        b_u = [b_u; 1+b_in-A_in(:,1:2)*x0'; 1-b_in+A_in(:,1:2)*x0'];

        % backward
        % A_x_ = [A_x_; A_x*Ctrl_m(:,1:2)];
        % b_x_ = [b_x_; b_x];
        % A_u = [A_u; A_in(:,1:2); -A_in(:,1:2)];
        % b_u = [b_u; u_max+b_in-A_in(:,3:4)*x0'; u_max-b_in+A_in(:,3:4)*x0'];
    end
    %%% Then add constraints such that x1, x2, and x3 satisfy the constraints
    %%% too
    A_f = [A_u; A_x_];
    b_f = [b_u; b_x_];

    V = Poly.hyp2vert(A_f, b_f);
    if size(V,1)>2
        V = Poly.conv(V);
        if earase_plots
            fR.XData = V(:,1);
            fR.YData = V(:,2);
        else
            patch(V(:,1),V(:,2),'g','facealpha',0.05);
        end
    else
        fR.XData = 0;
        fR.YData = 0;
    end

    %     drawnow;
    G.Nodes.x(ind_+1,:) = x_center;
    G.Nodes.F{ind_+1} = [A_f b_f];
    G.Nodes.B{ind_+1} = [A_b b_b];
    B_{ind_} = [A_b b_b];
    F_{ind_} = [A_f b_f];
    ind_ = ind_+1;
%     pause
end
G = rmnode(G,1);
drawnow;

%% Create Graph
ind = 1;
s = [];
t = [];
w = [];
for i = 1:G.numnodes
    F = F_{i};
    for j = 1:G.numnodes
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
G = G.addedge(s,t,w);
%% Plan a path
% start at the origin and end at a far eq. pt.
start_n = floor(pos_density/2)*vel_density + floor(vel_density/2)+1;
% end_n = (pos_density-1)*vel_density+ floor(vel_density/2)+1;
end_n = (pos_density-1)*floor(5*vel_density/6)+ floor(vel_density/2)+1-2;
scatter(G.Nodes.x(start_n,1), G.Nodes.x(start_n,2), 50, 'g', 'filled')
scatter(G.Nodes.x(end_n,1), G.Nodes.x(end_n,2), 50, 'b', 'filled')
path = shortestpath(G,start_n, end_n);
plot(G.Nodes.x(path,1),G.Nodes.x(path,2),'bo-','linewidth',2);

x_nodes = G.Nodes.x(path,:);
P = [];
V = [];
for i = 1:size(x_nodes,1)-1
    A_in = [G.Nodes.F{path(i)}(:,1:end-1); G.Nodes.B{path(i+1)}(:,1:end-1)];
    b_in = [G.Nodes.F{path(i)}(:,end); G.Nodes.B{path(i+1)}(:,end)];
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

%%
function M_ = projectOntoSemidefiniteCone(M)
% project onto semidefinite cone
[evec,eval] = eig(M);
M_ = eps*ones(size(M,1));
for i = 1:size(M,1)
    if eval(i,i) > 0
        M_ = M_ + evec(:,i)*evec(:,i)'*eval(i,i);
    end
end
end