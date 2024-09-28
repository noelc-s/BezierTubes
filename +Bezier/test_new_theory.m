% Parameters
u_max = 10;
horizon_N = 1;
steps = 1;
dt = 1;
A_x = [1 0; -1 0; 0 1; 0 -1];
b_x = [1; 1;1;1];

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
k_bar = 0;
X0 = x0;
XBAR = xbar;

e_0 = 0.05;

system.L.L_k = 1;
system.L.L_e = 0.0;
system.L.L_psi = 1;
system.L.L_pi = 1;
system.L.Lf = 1;
system.L.Lg = .0001;
% system.L.Lg = 1;
system.n = 2;
system.m = 1;
system.gamma = 2;
system.order = 3;

H = Bezier.H(system.order, horizon_N*dt);
Delta_vec = Bezier.Delta_vec(system.m, system.order, system.gamma);
H_vec = Bezier.H_vec(H, system.m, system.order, system.gamma, system.gamma-1);
D_vec = Delta_vec*H_vec;

% state constraints
Ax = system.L.L_pi*system.L.L_e*sqrt(diag(A_x*A_x'));
Ax = [zeros(size(Ax,1),1) Ax];
bx = b_x - system.L.L_pi*e_0*sqrt(diag(A_x*A_x'));

[Pi_x_1, Pi_x_2, delta_x] = Bezier.A_to_Pi(g_xbar, Ax, bx, system, H, xbar, f_xbar);

% input constraints
Au = [system.L.L_k*(1+system.L.L_psi) 1+system.L.L_k*system.L.L_e];
bu = u_max - (norm(k_bar, 2) + e_0);
[Pi_u_1, Pi_u_2, delta_u] = Bezier.A_to_Pi(g_xbar, Au, bu, system, H, xbar, f_xbar);

order = system.order;
gamma = system.gamma;
m = system.m;
j=1;

C_add = kron(A_x,ones(4*gamma*m*m,1));

Pi_input_constraints = [Pi_u_1 Pi_u_2];
Pi_state_constraints = [Pi_x_1 Pi_x_2];
Pi = [Pi_u_1 Pi_u_2; Pi_x_1+C_add Pi_x_2;];
h = [kron(bu,ones(4*gamma*m*m,1)) + Pi_input_constraints*[xbar(:,j); f_xbar(:,j)];
    kron(bx,ones(4*gamma*m*m,1)) + Pi_state_constraints*[xbar(:,j); f_xbar(:,j)]];
% Pi = [Pi_u_1 Pi_u_2; Pi_x_1 Pi_x_2;];
% h = [kron(delta_u,ones(4*gamma*m*m,1)) + Pi_input_constraints*[xbar(:,j); f_xbar(:,j)];
%     kron(delta_x,ones(4*gamma*m*m,1)) + Pi_state_constraints*[xbar(:,j); f_xbar(:,j)]];
% inds_where_Pi_is_zero = all(Pi_x_1==0,2);
% h = [1 + Pi_norm_constraints*[xbar(:,j); f_xbar(:,j)];
%     [(1-inds_where_Pi_is_zero).*(1-[Pi_x_1 Pi_x_2]*[xbar(:,j); f_xbar(:,j)]);
%     (1-inds_where_Pi_is_zero).*(1-[Pi_x_1 Pi_x_2]*[xbar(:,j); f_xbar(:,j)]);
%     (1-inds_where_Pi_is_zero).*(1-[Pi_x_1 Pi_x_2]*[xbar(:,j); f_xbar(:,j)]);
%     (1-inds_where_Pi_is_zero).*(1-[Pi_x_1 Pi_x_2]*[xbar(:,j); f_xbar(:,j)])]+...
%     [inds_where_Pi_is_zero*bx(1);
%     inds_where_Pi_is_zero*bx(2);
%     inds_where_Pi_is_zero*bx(3);
%     inds_where_Pi_is_zero*bx(4)]];

H_vec = Bezier.H_vec(H, m, order, gamma, gamma);
H_vec2 = Bezier.H_vec(H, m, order, gamma, gamma-1);
F = kron(eye(order+1),Pi)*H_vec;
G = kron(ones(order+1,1),h);

A = F;
b = G;

% Forward
Vert = cddmex('extreme',struct('A',[D_vec(1:2,:); A],'B',[X0;b],'lin',1:2));
Vert = Bezier.Poly.conv((D_vec(3:4,:)*Vert.V')');

if ~isempty(Vert)
    ind = convhull(Vert);
    Vert = Vert(ind,:);
end

clf
    axis off
    hold on;
    scatter(X0(1),X0(2),50,'filled','k')
    patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
    axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-.1 b_x(3)+.1]);
    patch(Vert(:,1),Vert(:,2),[1 0 0],'facealpha',0.03,'edgecolor',[0.7 0.1 0.1],'linewidth',2);
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'FontSize',17)
    set(gca,'linewidth',2)
    title('$u_{max} = 10$','interpreter','latex');
    axis square
    axis equal

%%
% tau = linspace(0,2*pi)';
% Circ = [cos(tau) sin(tau)];
% 
% num_interp = 20;
% for int = ((size(Vert,1)-1)*num_interp-1):-1:1
%     %     for i=1:size(Vert,1)
%     node = floor(int/num_interp)+1;
%     lambda = (node)*num_interp-int;
%     x1 = lambda/num_interp*Vert(node,:)'+(num_interp-lambda)/num_interp*Vert(node+1,:)';
%     %     x1 = [0;.9];
%     H = Bezier.H(order, dt);
%     D = Bezier.D(gamma,order, dt);
%     X = [x0 x1];
%     p = X(:)'/D;
%     P = [p; p*H];
%     
%     e = abs(system.L.L_e*p*H^2)+e_0;
%     
%     figure(1)
%     
%     clf
%     axis off
%     hold on;
%     scatter(X0(1),X0(2),50,'filled','k')
%     patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
%     axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-.1 b_x(3)+.1]);
%     patch(Vert(:,1),Vert(:,2),[1 0 0],'facealpha',0.03,'edgecolor',[0.7 0.1 0.1],'linewidth',2);
%     set(gca,'TickLabelInterpreter', 'latex');
%     set(gca,'FontSize',17)
%     set(gca,'linewidth',2)
%     title('$u_{max} = 10$','interpreter','latex');
%     axis square
%     axis equal
%     
%     
%     scatter(P(1,:),P(2,:),'filled')
%     Z = Bezier.Z(order, dt);
%     tau = linspace(0,1);
%     B = P*Z(tau);
%     hold on;
%     plot(B(1,:),B(2,:))
%     for k = 1:size(B,2)
%         plot(B(1,k)'+e*Z(tau(k))*Circ(:,1),...
%             B(2,k)'+e*Z(tau(k))*Circ(:,2),'k');
%     end
%     drawnow
% end
