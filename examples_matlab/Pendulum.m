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
u_max = 1;
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
K = [-1 -1];
w_overline = 0;

A_cl = [0 1; K];
P = lyap(A_cl',eye(2));
beta = 4*max(eig(P))^3/1^2;
state_buffer = b_x - sqrt(beta*w_overline^2)*sqrt(diag(A_x*inv(P)*A_x'));
e_bar = sqrt(beta*w_overline^2/min(eig(P)));

% Reference point 
x0 = [0; 0];
xbar = [0; 0];
f_xbar = f(xbar');
g_xbar = 1./g(xbar'); % This is g_inverse

[Q,Q_combined] = Bezier.Q(horizon_N, 3);

% Bezier Matrices
order = 3;
gamma = 2;
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order,dt);
D_nT = inv(D);
% Pi = Bezier.Pi(c,2,1);
Z = Bezier.Z(order, dt);
Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;

fig = figure(1);
clf;

plot_traj = false;
earase_plots = true;
plot_anything = false;

%%
% pos_density = 9;
% vel_density = 9;
pos_bounds = [-5 5];
vel_bounds = [-5 5];
% pos_bounds = [pi pi];
% vel_bounds = [0 0];
delta_pos = .5;
delta_vel = .5;

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
X_BAR_F = cell(numel(X),1);
X_BAR_B = cell(numel(X),1);
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
    x_bar_f = x_bar;
    f_xbar = f(x_bar')';
    g_xbar = g(x_bar')';
    [F,G] = Bezier.F_G(A_x, state_buffer, H,m, x_bar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
%     A_f = F(:,3:4);
%     b_f = G-F(:,1:2)*D_vec(1:2,1:2)*x_k;
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
    x_bar_b = x_bar;
    f_xbar = f(x_bar')';
    g_xbar = g(x_bar')';
    [F,G] = Bezier.F_G(A_x, state_buffer, H,m, x_bar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
%     A_b = F(:,1:2);
%     b_b = G-F(:,3:4)*D_vec(3:4,3:4)*x_k;
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

    if plot_anything
        drawnow;
    end
    Gr.Nodes.x(ind_+1,:) = x_k;
    Gr.Nodes.F{ind_+1} = [A_f b_f];
    Gr.Nodes.B{ind_+1} = [A_b b_b];
    X_K{ind_} = x_k;
    X_BAR_F{ind_} = x_bar_f;
    X_BAR_B{ind_} = x_bar_b;
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
                objective=[1 1];
                %         cddmex('extreme',struct('A',[D_vec(1:2,:); F],'B',[x_k;G],);
%                 IN= struct('obj',objective,'A',[D_vec;A_in],'B',[X_K{i};X_K{j};b_in],'lin',1:4);
                IN= struct('obj',objective,'A',[A_in],'B',[b_in]);
                OUT = cddmex('solve_lp',IN);
                if OUT.how==1
%                     Vert = Bezier.Poly.hyp2vert(A_in, b_in);
%                     if size(Vert,1)>2
%                         Vert = Bezier.Poly.conv(Vert);
                        s(ind) = i;
                        t(ind) = j;
                        if plot_anything && plot_traj
                            line([X_K{i}(1) X_K{j}(1)],[X_K{i}(2) X_K{j}(2)])
                        end
%                         w(ind) = norm(X_K{i} - mean(Vert))+norm(X_K{j}(1) - mean(Vert));
                        w(ind) = norm(X_K{i} - OUT.xopt)+norm(X_K{j}(1) - OUT.xopt);
%                         w(ind) = 1;
                        ind = ind+1;
%                     end
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
subplot(2,2,1)
scatter(Gr.Nodes.x(start_n,1), Gr.Nodes.x(start_n,2), 50, 'g', 'filled')
scatter(Gr.Nodes.x(end_n,1), Gr.Nodes.x(end_n,2), 50, 'b', 'filled')
path = shortestpath(Gr,start_n, end_n);
plot(Gr.Nodes.x(path,1),Gr.Nodes.x(path,2),'bo-','linewidth',2);

x_nodes = Gr.Nodes.x(path,:);
Pos = [];
Vel = [];
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
    Vel = [Vel C(2,:)];
    T = [T tau+T(end)];
    
    x0 = int_pt';
    x1 =  x_nodes(i+1,:)';
    X = [x0 x1];
    P = X(:)'/D;
    Xi = [P; P*H];
    tau = linspace(0,dt);
    C = Xi*Z(tau);
    Pos = [Pos C(1,:)];
    Vel = [Vel C(2,:)];
    T = [T tau+T(end)];
end
T = T(2:end);

%% MPC
addpath(genpath('~/mosek'))
addpath(genpath('~/repos/yalmip'))
disp('Setting up MPC problem');
yalmip('clear')

dt_MPC = dt/horizon_N;
N_MPC = 2*horizon_N;
state_stage_cost = 10;
input_stage_cost = 1;

%%% MPC size
state_dim = 2;
ctrl_pt_dim = order+1;
state_size_without_terminal = state_dim*(N_MPC);
state_size_with_terminal = state_dim*(N_MPC+1);
ctrl_pt_size_with_terminal = ctrl_pt_dim*N_MPC;
state_end_index=state_dim*N_MPC;
input_dim = 1;
input_size = input_dim*(N_MPC-1);
input_end_index=state_end_index+input_dim*(N_MPC-1);
total_dim = state_dim*N_MPC + input_dim*(N_MPC-1);
const_dim = 4*(4*m*gamma*m+size(A_x,1));

% Pre-allocate program
X_MPC = [];
T_MPC = [];

% Bezier stuff
H = Bezier.H(order, dt_MPC);
D = Bezier.D(gamma,order,dt_MPC);
Z = Bezier.Z(order, dt_MPC);
Delta_vec = Bezier.Delta_vec(m, order, gamma);
H_vec = Bezier.H_vec(H, m, order, gamma, gamma-1);
D_vec = Delta_vec*H_vec;
Q = {eye(4)};

% Optimizer variables
ctrl_pts = sdpvar(ctrl_pt_dim,N_MPC);
states = sdpvar(state_dim,N_MPC+1);
states_ref_ = sdpvar(state_dim,N_MPC+1);
% state_error = sdpvar(state_dim,1);
% inputs = sdpvar(input_dim,N_MPC-1);
x0_ = sdpvar(state_dim,1);
xf_ = sdpvar(state_dim,1);
F_ = sdpvar(const_dim, ctrl_pt_size_with_terminal);
G_ = sdpvar(const_dim, N_MPC);

Constraints = [];
Constraints = [Constraints states(:,1) == x0_];
% Constraints = [Constraints states(:,end) == xf_];
for i = 1:N_MPC
    Constraints = [Constraints F_(:,(i-1)*ctrl_pt_dim+1:i*ctrl_pt_dim)*ctrl_pts(:,i)...
        <= G_(:,i)];
    Constraints = [Constraints D_vec*ctrl_pts(:,i) == [states(:,i); states(:,i+1)]];
end

Q_cost = state_stage_cost*eye(state_dim*(N_MPC+1));
R = input_stage_cost*eye(input_dim*(N_MPC-1));

state_error = states - states_ref_;

Objective = 1/2*(state_error(:)'*Q_cost*state_error(:)); % change if you want different state weighting
P = optimizer(Constraints,Objective,sdpsettings('solver','mosek','verbose',0),...
    {x0_,xf_,states_ref_,F_,G_},{states,ctrl_pts});

x0 = x_nodes(1,:)';
F = zeros(const_dim, ctrl_pt_size_with_terminal);
G = zeros(const_dim, N_MPC);

x_bar = [];
for i = 1:size(x_nodes,1)-1
    x_bar = [x_bar X_BAR_F{path(i)}];
    x_bar = [x_bar X_BAR_B{path(i+1)}];
end

for i = 1:size(x_bar,2)-N_MPC
    for j = i:(N_MPC+i-1)
        xbar = x_bar(:,j);
        f_xbar = f(xbar')';
        g_xbar = g(xbar')';
        [F_tmp,G_tmp] = Bezier.F_G(A_x, state_buffer, H,m, xbar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
        F(:,((j-i+1)-1)*ctrl_pt_dim+1:(j-i+1)*ctrl_pt_dim) = F_tmp;
        G(:,(j-i+1)) = G_tmp;
    end
    
    X = [Pos; Vel]';
    [~,ind_T] = unique(T);
    states_ref = interp1(T(ind_T),X(ind_T,:),dt_MPC*(i-1) + (0:dt_MPC:N_MPC*dt_MPC))';

% X_ = x_nodes;
% T_ = 0:dt:dt*(size(x_nodes,1)-1);
% states_ref = interp1(T_,X_,dt_MPC*(i-1) + (0:dt_MPC:N_MPC*dt_MPC))';

    xf = [1;1]; % random for now
    
    [sol, diagnostics,d1,d2,d3,d4] = P({x0,xf,states_ref,F,G});
    X_MPC = [X_MPC; sol{1}(:,1)'];
    T_MPC = [T_MPC dt_MPC*(i-1)];
    x0 = sol{1}(:,2);
end


% for i = 1:size(x_nodes,1)-2
%     for ind = 0:N_MPC-1
%         if (ind+1)>=N_MPC/2
%             set = 2;
%         else
%             set = 1;
%         end
%     indeces = ind+(1:N_MPC);
%     for j = 1:length(indeces)
%         index = mod(indeces(j),N_MPC/2)+1;
%         if index == 1
%            set = set+1; 
%         end
%         if set == 1
%             x_bar(:,j) = X_BAR_F{path(i)}(:,index);
%         elseif set == 2
%             x_bar(:,j) = X_BAR_B{path(i+1)}(:,index);
%         elseif set == 3
%             x_bar(:,j) = X_BAR_F{path(i+1)}(:,index);
%         end
%         
%         f_xbar = f(x_bar(:,j)')';
%         g_xbar = g(x_bar(:,j)')';
%         [F_tmp,G_tmp] = Bezier.F_G(A_x, b_x, H,m, x_bar(:,j), f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
%         F(:,(j-1)*ctrl_pt_dim+1:j*ctrl_pt_dim) = F_tmp;
%         G(:,j) = G_tmp;
%     end
%     
% %     for j = 1:N_MPC/2
% %         x_bar = X_BAR_F{path(i)}(:,j);
% %         f_xbar = f(x_bar')';
% %         g_xbar = g(x_bar')';
% %         [F_tmp,G_tmp] = Bezier.F_G(A_x, b_x, H,m, x_bar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
% %         F(:,(j-1)*ctrl_pt_dim+1:j*ctrl_pt_dim) = F_tmp;
% %         G(:,j) = G_tmp;
% %     end
% %     for j = 1:N_MPC/2
% %         x_bar = X_BAR_B{path(i+1)}(:,j);
% %         f_xbar = f(x_bar')';
% %         g_xbar = g(x_bar')';
% %         [F_tmp,G_tmp] = Bezier.F_G(A_x, b_x, H,m, x_bar, f_xbar, g_xbar, 2,Q,Lg,Lf,e_bar,K,u_max);
% %         F(:,(j-1)*ctrl_pt_dim+N_MPC/2*ctrl_pt_dim+1:j*ctrl_pt_dim+N_MPC/2*ctrl_pt_dim) = F_tmp;
% %         G(:,j+N_MPC/2) = G_tmp;
% %     end
% %         states_ref = repmat(x_nodes(i+1,:)',1,N_MPC+1);
%     %     states_ref = [x0 X_BAR_F{path(i)} X_BAR_B{path(i+1)}];
%     X = [Pos; Vel]';
%     [~,ind_T] = unique(T);
%     states_ref = interp1(T(ind_T),X(ind_T,:),dt_MPC*ind+(2*dt*(i-1):dt_MPC:2*dt*i))';
%     %     x0 = x_nodes(i,:)';
%     xf = [1;1]; % random for now
%     
%     [sol, diagnostics,d1,d2,d3,d4] = P({x0,xf,states_ref,F,G});
%     X_MPC = [X_MPC; sol{1}(:,1)'];
%     T_MPC = [T_MPC dt_MPC*ind+2*dt*(i-1)];
%     x0 = sol{1}(:,2);
%     end
% end

% X_MPC = [X_MPC; x0'];
clear set


%% Plot things

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

pause
subplot(2,2,[3 4]); hold on;
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)
plot(T,Pos,'linewidth',2);
plot(T,Vel,'linewidth',2);
plot(T_MPC, X_MPC,'linewidth',2)

% Animate
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