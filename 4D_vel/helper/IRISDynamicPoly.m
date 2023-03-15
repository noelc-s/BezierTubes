%% Algorithm 1 for 2d


if ~exist('obstacle_number')
    obstacle_number = 10;
end

O = {};

max_vel = 10;

A_x = [1 0 0 0; -1 0 0 0; 0 1 0 0; 0 -1 0 0;... % pos
    0 0 1 0; 0 0 0 1; 0 0 -1 0; 0 0 0 -1];      % vel
b_x = [2; 2; 2; 2;...
    max_vel; max_vel; max_vel; max_vel];

% box_r = 0.1;
% for i = 1:obstacle_number
%     center  = rand(2,1)*3-1.5;
%     %     center = [0 i/2-2.5]';
%     O{i} = box_r*[-1 1; -1 -1; 1 -1; 1 1]'+center;
% end

O{1} = 1*[-1 1 max_vel max_vel; -1 -1 max_vel max_vel; 1 -1 max_vel max_vel; 1 1 max_vel max_vel;...
    -1 1 -max_vel max_vel; -1 -1 -max_vel max_vel; 1 -1 -max_vel max_vel; 1 1 -max_vel max_vel;...
    -1 1 max_vel -max_vel; -1 -1 max_vel -max_vel; 1 -1 max_vel -max_vel; 1 1 max_vel -max_vel;...
    -1 1 -max_vel -max_vel; -1 -1 -max_vel -max_vel; 1 -1 -max_vel -max_vel; 1 1 -max_vel -max_vel]'+[0; 0; 0; 0];

% O{2} = 1*[-1 1 max_vel max_vel; -1 -1 max_vel max_vel; 1 -1 max_vel max_vel; 1 1 max_vel max_vel;...
%     -1 1 -max_vel max_vel; -1 -1 -max_vel max_vel; 1 -1 -max_vel max_vel; 1 1 -max_vel max_vel;...
%     -1 1 max_vel -max_vel; -1 -1 max_vel -max_vel; 1 -1 max_vel -max_vel; 1 1 max_vel -max_vel;...
%     -1 1 -max_vel -max_vel; -1 -1 -max_vel -max_vel; 1 -1 -max_vel -max_vel; 1 1 -max_vel -max_vel]'+[0; .8; 0; 0];

for i = 1:size(O,2)
    [A_iris,b_iris] = vert2lcon(O{i}');
    A_O{i} = A_iris;
    b_O{i} = b_iris;
end

if ~exist('density')
    density = 5;
end
if ~exist('overlap')
    overlap = false;
end
[X0, Y0] = meshgrid(linspace(-1.95,1.95,density));

ind = 1;

hold on;
line([-2 -2], [2 -2],'color','k')
line([2 2], [2 -2],'color','k')
line([-2 2], [2 2],'color','k')
line([-2 2], [-2 -2],'color','k')

for i = 1:size(O,2)
    patch(O{i}(1,:),O{i}(2,:),'r','facealpha',0.05)
end
axis equal

Polytopes = {};
PolyCenter = {};
X0 = X0(:);
Y0 = Y0(:);
X0(end+1) = IC(1);
Y0(end+1) = IC(2);
X0(end+1) = EC(1);
Y0(end+1) = EC(2);
for t = 1:numel(X0)
    q0 = [X0(t) Y0(t) 0 0]';
    cont = true;

    if t < numel(X0)-1
        for j = 1:size(O,2)
            if A_O{j}*q0 <= b_O{j}+0.01 % add some tolerance so if yo usample near boundary it doesn't stall
                cont = false;
                break;
            end
        end
        if overlap
            for p = Polytopes
                if p{1}(:,1:end-1)*q0 <= p{1}(:,end)
                    cont = false;
                    break;
                end
            end
        end
        if ~cont
            continue
        end
    end

    stopping_tol = 1e-5;

    eps = 1e-3;
    C_i = 10*eps*eye(4);
    C_ip1 = eps*eye(4);
    d_ip1 = q0;
    i = 0;

    A_ip1 = [];
    b_ip1 = [];

    [A_dyn, b_dyn] = Poly.dynamicTube2D(H, A, B,A_x, b_x,u_max,D_nT);
    
%     [A_dyn_R, b_dyn_R] = Poly.backwardReachable2D(H, A, B,A_x, b_x,u_max,D_nT,[q0]');

%     A_dyn = [A_dyn_F; A_dyn_R];
%     b_dyn = [b_dyn_F; b_dyn_R];


    % Should not do forward and backward reachable, this is artificially
    % restrictive.
    b_dyn_F = b_dyn - A_dyn(:,1:4)*q0; 
    A_dyn_F = A_dyn(:,5:8);

    b_dyn_R = b_dyn - A_dyn(:,5:8)*q0; 
    A_dyn_R = A_dyn(:,1:4);

    A_dyn = [A_dyn_F; A_dyn_R];
    b_dyn = [b_dyn_F; b_dyn_R];

    A_bounds = A_dyn;
    b_bounds = b_dyn;

    while(abs((det(C_ip1) - det(C_i))/det(C_i)) > stopping_tol)

        C_i = C_ip1;
        d_i = d_ip1;
        A_i = A_ip1;
        b_i = b_ip1;

        [A_ip1, b_ip1] = IRIS.separatingHyperplanes_cell(C_i, d_i, O);

        A_ip1 = [A_ip1; A_bounds];
        b_ip1 = [b_ip1; b_bounds];

        [C_ip1, d_ip1] = IRIS.inscribedEllipsoid(A_ip1, b_ip1);
        i = i+1;
    end

    A_iris = A_ip1;
    b_iris = b_ip1;
    C = C_ip1;
    d = d_ip1;

    %     if (A_dyn*q0 <= b_dyn) % don't even consider the point if it is not in the iris polytope
    A_dyn = A_iris;
    b_dyn = b_iris;

    add = true;
    for i = 1:length(Polytopes)
        if size(Polytopes{i},1) ~= 0
            nonzero_ind = ~(sum(Polytopes{i}(:,1:2)==0,2)==2);
            if size([A_dyn b_dyn],1) == size(Polytopes{i},1) & norm([A_dyn(nonzero_ind,:)] - Polytopes{i}(nonzero_ind,1:end-1))<=1e-2 & norm([A_dyn(nonzero_ind,:)] - Polytopes{i}(nonzero_ind,end))<=1e-2
                add = false;
                break
            end
        end
    end

    scatter(q0(1), q0(2),30,[50 190 50]/255,'filled')

    if add
        Polytopes{ind} = [A_dyn b_dyn];
        V = lcon2vert(A_dyn,b_dyn);
        PolyCenter{ind} = mean(V);
        if ~isempty(V)
            c_h = convhull(V(:,1:2));
            c_h = V(c_h,:);
            patch(c_h(:,1),c_h(:,2),'b','facealpha',0.01)
            scatter(PolyCenter{ind}(1), PolyCenter{ind}(2),100,'k','filled')
            text(PolyCenter{ind}(1), PolyCenter{ind}(2)+0.1, string(ind));
            ind = ind+1;
            if ~overlap
                O{size(O,2)+1} = V'; % add grown polytope to obstacle list
                [A_dyn,b_dyn] = vert2lcon(O{size(O,2)}');
                A_O{size(O,2)} = A_dyn;
                b_O{size(O,2)} = b_dyn;
            end
        end
    end
    drawnow
end
