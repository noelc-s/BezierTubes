addpath('..')
plot_traj = false;
earase_plots = false;

pos_density = 11;
vel_density = 11;
u_max = .2;
dt = 1;

A = [0 1; -1 0];
B = [0; 1];
S = LinearSystem(A,B);

A_x = [1 0; -1 0; 0 1; 0 -1];
b_x = [1; 1; 1; 1];

order = 2*size(A,1)-1;
[H, D_nT] = Poly.getBezMatrices(order, dt);

f = figure(1);
clf;

% subplot(1,2,1)
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

axis([-1.1 1.1 -1.1 1.1])
axis equal
V_x = Poly.conv(Poly.hyp2vert(A_x,b_x));
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
tic
for i = 1:numel(X)

    x_center = [X(i) Y(i)];

    %%% Forward Reachable
    x1 = x_center;
    [A_b, b_b] = Poly.backwardReachable(H, A, B,A_x, b_x,u_max,D_nT,x1);
    V = Poly.hyp2vert(A_b, b_b);
    if size(V,1)>2
        V = Poly.conv(V);
        if earase_plots
            bR.XData = V(:,1);
            bR.YData = V(:,2);
        else
            patch(V(:,1),V(:,2),'b','facealpha',0.05);
        end
        if plot_traj
            delete(t_bR);
            t_bR = Poly.plotReachableTraj(V,'b',[],x1,dt,H);
        end
    else
        bR.XData = 0;
        bR.YData = 0;
    end

    %%% Reverse Reachable
    x0 =x_center;
    [A_f, b_f] = Poly.forwardReachable(H, A, B,A_x, b_x,u_max,D_nT,x0);
    V = Poly.hyp2vert(A_f, b_f);
    if size(V,1)>2
        V = Poly.conv(V);
        if earase_plots
            fR.XData = V(:,1);
            fR.YData = V(:,2);
        else
            patch(V(:,1),V(:,2),'g','facealpha',0.05);
        end
        if plot_traj
            delete(t_fR);
            t_fR = Poly.plotReachableTraj(V,'g',x0,[],dt,H);
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
    i/G.numnodes
    for j = 1:G.numnodes
        B = B_{j};
        A_in = [F(:,1:end-1); B(:,1:end-1)];
        b_in = [F(:,end); B(:,end)];
        quadprog(
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
end_n = (pos_density-1)*vel_density+ floor(vel_density/2)+1;

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
    [~, ~, pos, vel] = Poly.plotTraj('r',x_nodes(i,:),int_pt,dt,H);
    plot(pos, vel,'linewidth',1,'color','r');
    P = [P pos(1:end-1)];
    V = [V vel(1:end-1)];
    [~, ~, pos, vel] = Poly.plotTraj('r',int_pt,x_nodes(i+1,:),dt,H);
    plot(pos, vel,'linewidth',1,'color','r');
    P = [P pos(1:end-1)];
    V = [V vel(1:end-1)];
end
drawnow
% end
toc
% subplot(1,2,2); hold on;
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)
plot(P,'linewidth',2);
plot(V,'linewidth',2);