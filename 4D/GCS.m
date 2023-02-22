%% Search graph for path
addpath('helper')
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

plot([path(:,1); EC(1)],[path(:,2); EC(2)],'ko--')
drawnow;