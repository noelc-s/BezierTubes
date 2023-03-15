%% Search graph for path
addpath('helper')
addpath(genpath('~/mosek'))
addpath('..')
clear;clf;
axis equal
axis off
tic

density = 50;
obstacle_number = 5;
dt = 5;
dt_short = 0.25;
N = 2*ceil(dt/dt_short);
u_max = 1;
overlap = false;

A = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
B = [0 0; 0 0; 1 0; 0 1];
order = 3;
[H, D_nT] = Poly.getBezMatrices(order, dt);

f = figure(1);
tg = uitabgroup(f);
t1 = uitab(tg, 'Title', 'RRT');
a1 = axes('Parent', t1);

IC = [-1.9 -1.9];
EC = [1.9 1.9];

IRISDynamicPoly;
BuildIntersectionGraphFromPoly;
convertDualGraphtoPython;

if overlap
    buffer = 0.0;
    PC = cell2mat(PolyCenter');
    [~,SI] = min(vecnorm(IC'-PC'));
    [~,EI] = min(vecnorm(EC'-PC'));
else
    buffer = 0.01;
    for i=1:size(Polytopes,2)
        if Polytopes{i}(:,1:end-1)*IC' <= Polytopes{i}(:,end)
            SI = i;
        end
        if Polytopes{i}(:,1:end-1)*EC' <= Polytopes{i}(:,end)
            EI = i;
        end
    end
end

indS = 1;
indE = 1;

scatter(IC(1), IC(2),50,'r','filled')
scatter(EC(1), EC(2),50,'g','filled')

%%
[status,cmdout] = system('source /home/amber/env/bin/activate; python3 helper/GraphConvexSets.py');
if status ~= 0
    error('Python errored')
end
gcs_optimized_path;
ord = [];
node = 0;
poly_num = SI;
for i = 1:size(edgeTraversal,1)
    ind = find(edgeTraversal(:,1)==node);
    ord = [ord ind];
    node = edgeTraversal(ind,2);
end
edgeTraversal = edgeTraversal(ord,1);
path = path(ord,1:2);

edgeTraversal = edgeTraversal-1;
edgeTraversal = edgeTraversal(2:end);

Ptopes = [G.Edges.EndNodes(edgeTraversal(1),1); G.Edges.EndNodes(edgeTraversal(1),2)];
node = G.Edges.EndNodes(edgeTraversal(1),2);
inds = 2:length(edgeTraversal);
for i = 1:length(edgeTraversal)
    row = find(sum(G.Edges.EndNodes(edgeTraversal(inds),:)==node,2));
    col = find(G.Edges.EndNodes(edgeTraversal(inds(row)),:)==node);
    node = G.Edges.EndNodes(edgeTraversal(inds(row)),mod(col,2)+1);
    tmp_inds = true(size(inds,2),1);
    tmp_inds(row) = 0;
    inds = inds(tmp_inds);
    Ptopes = [Ptopes; node];
end

path = [path; EC];

plot(path(:,1),path(:,2),'ko--')
drawnow;

%%

P = [];
Time = 0;
for i = 1:size(path,1)-1
[~, tau, pos_1, vel_1] = Poly.plotTraj('r',[path(i,1) 0],[path(i+1,1) 0],dt,H);
[~, ~, pos_2, vel_2] = Poly.plotTraj('r',[path(i,2) 0],[path(i+1,2) 0],dt,H);
P = [P;pos_1' pos_2'];
tau = tau*dt;
Time = [Time; tau'+Time(end)];
end
Time = Time(2:end);
[~,ind] = unique(Time);
Time = Time(ind);
P = P(ind,:);

new_time = linspace(Time(1),Time(end),1000);
P = interp1(Time,P,new_time);
Time = new_time;


s = scatter(0,0,50,'b','filled');
tic
i = 1;
while i < length(P)
    i = find(Time>toc,1);
    s.XData = P(i,1);
    s.YData = P(i,2);
    drawnow
end
%%

vel = [0 0];
PATH = [];
for i = 1:size(path,1)-1
x0 = path(i,:);
x1 = path(i+1,:);
opts = optimset('Display','on');
A_x = [Polytopes{Ptopes(i)}(:,1:end-1) zeros(size(Polytopes{Ptopes(i)},1),2)];
b_x = Polytopes{Ptopes(i)}(:,end);
% V = Poly.hyp2vert(A_x,b_x);
% V = [[V(:,1:2) .1*ones(size(V,1),2)];[V(:,1:2) -.1*ones(size(V,1),2)]];
% [A_x, b_x] =  Poly.vert2hyp(V);
% A_x = A_x/norm(b_x);
% b_x = b_x/norm(b_x)+0.01;
% [t,FVAL,EXITFLAG] = fmincon(@(t) t(1), [1 0 0],[],[],[],[],[],[],@(t) con(t, order,A,B,u_max,A_x,b_x,x0,x1),opts)
[t,FVAL,EXITFLAG] = fmincon(@(t) t(1), [1],[],[],[],[],[],[],@(t) con(t, order,A,B,u_max,A_x,b_x,x0,x1),opts)
T(i) = t(1);
vel = t(2:end);
% PATH = [PATH; x0 x1 vel];
end

P = [];
Time = 0;
for i = 1:size(path,1)-1
[~, tau, pos_1, vel_1] = Poly.plotTraj('r',[path(i,1) 0],[path(i+1,1) 0],T(i),H);
[~, ~, pos_2, vel_2] = Poly.plotTraj('r',[path(i,2) 0],[path(i+1,2) 0],T(i),H);
P = [P;pos_1' pos_2'];
tau = tau*T(i);
Time = [Time;tau'+Time(end)];
end

Time = Time(2:end);
[~,ind] = unique(Time);
Time = Time(ind);
P = P(ind,:);

new_time = linspace(Time(1),Time(end),1000);
P = interp1(Time,P,new_time);
Time = new_time;

s = scatter(0,0,50,'b','filled');
tic
i = 1;
while i < length(P)
    i = find(Time>toc,1);
    s.XData = P(i,1);
    s.YData = P(i,2);
    drawnow
end



function [c, ceq] = con(t, order,A,B,u_max,A_x,b_x,x0,x1)
[H, D_nT] = Poly.getBezMatrices(order, t(1));
H_0 = H^0; 
H_1 = H^1;
H_2 = H^2;
A_u = [];
b_u = [];
A_x_ = [];
b_x_ = [];
for m = 1:4
    I_m = zeros(2,8);
    I_m(1,(m-1)*2+1) = 1;
    I_m(2,(m-1)*2+2) = 1;
    A_tmp = [I_m*kron(H_0,eye(2))'; I_m*kron(H_1,eye(2))'];
    A_lin = ([1 0 0 0; 0 1 0 0]*A*B)\([1 0 0 0; 0 1 0 0]*A*A*A_tmp - I_m*kron(H_2,eye(2))');
    A_in = A_lin*kron(D_nT,eye(2));

    A_u = [A_u; A_in; -A_in];
    b_u = [b_u; [u_max; u_max; u_max; u_max]]; % two control inputs

    A_x_ = [A_x_; A_x*A_tmp*kron(D_nT,eye(2))];
    b_x_ = [b_x_; b_x];
end
A_ = [A_u; A_x_];
b_ = [b_u; b_x_];
c = A_*[x0(1) x0(2) 0 0 x1(1) x1(2) 0 0]' - b_;
% c = A_*[x0(1) x0(2) x0(3) x0(4) x1(1) x1(2) t(2) t(3)]' - b_;
ceq = [];
end