
plot_traj = false;
earase_plots = false;

u_max = .5;
dt = 1;
N = 17;

A = [0 1; -1 0];
B = [0; 1];
S = LinearSystem(A,B);

A_x = [1 0; -1 0; 0 1; 0 -1];
b_x = [4; 4; 4; 4];

order = 2*size(A,1)-1;
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

axis([-4.4 4.1 -4.1 4.1])
axis equal
V_x = lcon2vert(A_x,b_x);
ind_ux = convhull(V_x);
V_x = V_x(ind_ux,:);
patch(V_x(:,1),V_x(:,2),'k','facealpha',0.05,'linewidth',2);

set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)

% [X,Y] = meshgrid(linspace(-b_x(2)+0.05, b_x(1)-0.05,pos_density), linspace(-b_x(4)+0.05, b_x(3)-0.05, vel_density));
% 
% G = digraph;
% G = G.addnode(1);
% G.Nodes.x = [0 0];
% G.Nodes.F = {0};
% G.Nodes.B = {0};

% ind_ = 1;
% G = G.addnode(numel(X));
% B_ = cell(numel(X),1);
% F_ = cell(numel(X),1);

x0 = [0 0];
xT = [3.14 0];
[A_in, b_in] = Poly.dynamicTube(H, A, B,A_x, b_x,u_max,D_nT);


A_big = zeros(size(A_in,1)*(N-1),2*N);
b_big = [];
for i = 1:N-1
    A_big(((i-1)*size(A_in,1)+1):((i)*size(A_in,1)), ((i-1)*2+1):((i-1)*2+4)) = A_in;
    b_big = [b_big; b_in];
end

[x,fval,exitflag] = quadprog(zeros(2*N), zeros(2*N,1), A_big, b_big, [eye(2) zeros(2,2*(N-1)); zeros(2,2*(N-1)) eye(2)], [x0'; xT'])

x = [x(1:2:end) x(2:2:end)];
scatter(x(:,1),x(:,2),'filled')

P = [];
V = [];

for i = 1:size(x,1)-1
%     A_in = [G.Nodes.F{path(i)}(:,1:end-1); G.Nodes.B{path(i+1)}(:,1:end-1)];
%     b_in = [G.Nodes.F{path(i)}(:,end); G.Nodes.B{path(i+1)}(:,end)];
%     Vert = Poly.conv(Poly.hyp2vert(A_in, b_in));
%     int_pt = mean(Vert);
    [~, ~,pos, vel] = Poly.plotTraj('r',x(i,:),x(i+1,:),dt,H);
    P = [P pos(1:end-1)];
    V = [V vel(1:end-1)];
    [A_in, b_in] = Poly.forwardReachable(H, A, B,A_x, b_x,u_max,D_nT,x(i,:));
    Vert = Poly.hyp2vert(A_in, b_in);
    Vert = Poly.conv(Vert);
    patch(Vert(:,1),Vert(:,2),'g','facealpha',0.05,'linewidth',.5);
    if i < size(x,1)-1
    [A_inb, b_inb] = Poly.backwardReachable(H, A, B,A_x, b_x,u_max,D_nT,x(i+2,:));
    Vert = Poly.hyp2vert(A_inb, b_inb);
    Vert = Poly.conv(Vert);
    patch(Vert(:,1),Vert(:,2),'b','facealpha',0.05,'linewidth',.5);
    Vert = Poly.hyp2vert([A_in; A_inb], [b_in; b_inb]);
    if ~isempty(Vert)
    Vert = Poly.conv(Vert);
    patch(Vert(:,1),Vert(:,2),'r','facealpha',0.5,'linewidth',2);
    end
    end
end
subplot(1,2,2); hold on;
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'FontSize',17)
set(gca,'linewidth',2)
plot(linspace(0,N*dt,length(P)),P,'linewidth',2);
plot(linspace(0,N*dt,length(P)),V,'linewidth',2);