%% Description
% A script to test compare three things:
% - N step reachable over dt
% - 1 step reachable over N*dt
% - 1 step reachable with N step B-spline over N*dt

% Parameters
u_max =5;
horizon_N = 10;
dt = .1;
A_x = [0 1; 0 -1; 1 0; -1 0];
b_x = 4*[1; 1; 1; 1];

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
x0 = [0; 1];

xbar = x0;
f_xbar = f(xbar');
g_xbar = g(xbar');
X0 = x0;
XBAR = xbar;

Lf = 1;
Lg = 1;
e_bar = 0;
K = [-1 -1];
[M, N, Gamma, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);
clf

%% 1 step reachable over N*dt
% Bezier Matrices for N*dt
H = Bezier.H(3, horizon_N*dt);
D = Bezier.D(2,3,horizon_N*dt);
D_nT = inv(D');
Pi = Bezier.Pi(c,2,1);
Z = Bezier.Z(3, horizon_N*dt);

color = 'g';
tic
[A_in, b_in] = Bezier.F_G(A_x, b_x, Pi, H, xbar, f_xbar, 2);
A = A_in; 
b = b_in;

IC = x0;

x0 = IC;


% Forward
Vert = cddmex('extreme',struct('A',[D(1:2,:); A],'B',[x0;b],'lin',1:2));
Vert = Vert.V*D(3:4,:)';
toc

% Vert = Bezier.Poly.hyp2vert(A,b);

if ~isempty(Vert)
ind = convhull(Vert);
Vert = Vert(ind,:);
end
%     patch(Vert(:,1),Vert(:,2),color,'facealpha',0.1);
% end
Vert1 = Vert;

%% 1 step B-spline reachable over dt, different x_bar
% Bezier Matrices for N*dt
H = Bezier.H(3, horizon_N*dt);
D = Bezier.D(2,3,horizon_N*dt);
D_nT = inv(D');
Pi = Bezier.Pi(c,2,1);
Z = Bezier.Z(3, horizon_N*dt);

x1 = X0;
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

xbar = x_bar;
f_xbar = f(x_bar')';

Q = Bezier.Q(horizon_N, 3);
tic
[A_in, b_in] = Bezier.F_G(A_x, b_x, Pi, H, xbar, f_xbar, 2,Q);
A = A_in; 
b = b_in;

% Forward
Vert = cddmex('extreme',struct('A',[D(1:2,:); A],'B',[X0;b],'lin',1:2));
Vert = Vert.V*D(3:4,:)';
toc
% 
% H_0 = H^0;
% H_1 = H^1;
% H_2 = H^2;
% A_x_ = [];
% b_x_ = [];
% A_u_ = [];
% b_u_ = [];
% 
% tic
% for m = 1:4
%     for k = 1:length(Q)
%         I_m = zeros(1,4);
%         I_m(m) = 1;
%         A_tmp = [I_m*Q{k}*H_0'; I_m*Q{k}*H_1'];
%         
%         x_bar = xbar(:,k);
%         f_xbar = f(x_bar');
%         g_xbar = g(x_bar');
%         [M, N, Gamma, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);
%         Pi = Bezier.Pi(c,2,1);
%         
%         % input constaints
%         
%         A_u_ = [A_u_; Pi*[eye(2)*A_tmp*D_nT;...
%             I_m*Q{k}*H_2'*D_nT]];
%         b_u_ = [b_u_; 1+Pi*[x_bar; f_xbar]]; % has to match numerator of c.
%         
%         
%         % and state constraints
%         % If this were over entire horizon, it would need different H
%         A_x_ = [A_x_; A_x*A_tmp*D_nT];
%         b_x_ = [b_x_; b_x];
%     end
% end
% 
% A_in = [A_u_; A_x_];
% b_in = [b_u_; b_x_];
% 
% IC = X0;

% A = A_in(:,3:4);
% b = b_in - A_in(:,1:2)*IC;
% toc

% subplot(3,2,6)
% hold on;
% scatter(xbar(1,:),xbar(2,:))
% patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
% axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-0.1 b_x(3)+0.1]);
% Vert = Bezier.Poly.hyp2vert(A,b);
if ~isempty(Vert)
ind = convhull(Vert);
Vert = Vert(ind,:);
end
% patch(Vert(:,1),Vert(:,2),color,'facealpha',0.1);
Vert2 = Vert;

%% N step reachable over dt, same x_bar
% Bezier Matrices for dt

% Bezier Matrices for dt

H = Bezier.H(3, dt);
D = Bezier.D(2,3,dt);
D_nT = inv(D');
Pi = Bezier.Pi(c,2,1);
Z = Bezier.Z(3, dt);
H_0 = H^0;
H_1 = H^1;
H_2 = H^2;

% subplot(3,2,4)
% hold on;

% patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
% axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-0.1 b_x(3)+0.1]);

Vert_ = X0';

        xbar = IC;
        f_xbar = f(xbar');
        g_xbar = g(xbar');
        [M, N, Gamma, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);
        Pi = Bezier.Pi(c,2,1);

tic

for k = 1:horizon_N
    Vert_tmp = [];
    for i = 1:size(Vert_,1)
        IC = Vert_(i,:)';
        
%         
%         A_x_ = [];
%         b_x_ = [];
%         A_u_ = [];
%         b_u_ = [];
%         for m = 1:4
%             I_m = zeros(1,4);
%             I_m(m) = 1;
%             A_tmp = [I_m*H_0'; I_m*H_1'];
%             
%             % input constaints
%             A_u_ = [A_u_; Pi*[eye(2)*A_tmp*D_nT;...
%                 I_m*H_2'*D_nT]];
%             b_u_ = [b_u_; 1+Pi*[xbar; f_xbar]]; % has to match numerator of c.
%             
%             % and state constraints
%             A_x_ = [A_x_; A_x*A_tmp*D_nT];
%             b_x_ = [b_x_; b_x];
%         end
%         A_in = [A_u_; A_x_];
%         b_in = [b_u_; b_x_];
%         
%         A = A_in(:,3:4);
%         b = b_in - A_in(:,1:2)*IC;

        [A_in, b_in] = Bezier.F_G(A_x, b_x, Pi, H, xbar, f_xbar, 2);
        A = A_in; 
        b = b_in;

        
        x0 = IC;
        color = 'g';
                Vert = cddmex('extreme',struct('A',[D(1:2,:); A],'B',[x0;b],'lin',1:2));
        Vert = uniquetol(Vert.V*D(3:4,:)',1e-3,'ByRows',1);
        
%         try
%             Vert = uniquetol(Bezier.Poly.hyp2vert(A,b),1e-3,'ByRows',1);
%         catch E
%             break
%         end
        if ~isempty(Vert)
            if size(Vert,1) > 2
                ind = convhull(Vert);
                
                Vert = Vert(ind,:);
            end
            Vert_tmp = [Vert_tmp; Vert(1:end-1,:)];
        end
        patch(Vert(:,1),Vert(:,2),color,'facealpha',0.1);
    end
    Vert_ = Vert_tmp;
    ind = convhull(Vert_);
    Vert_ = Vert_(ind,:);
        
end
toc
Vert3 = Vert_;

%% N step reachable over dt, different x_bar
% Bezier Matrices for dt

H = Bezier.H(3, dt);
D = Bezier.D(2,3,dt);
D_nT = inv(D');
Pi = Bezier.Pi(c,2,1);
Z = Bezier.Z(3, dt);
H_0 = H^0;
H_1 = H^1;
H_2 = H^2;

% subplot(3,2,4)
% hold on;

% patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
% axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-0.1 b_x(3)+0.1]);

Vert_ = X0';
tic
for k = 1:horizon_N
    Vert_tmp = [];
    for i = 1:size(Vert_,1)
        IC = Vert_(i,:)';
        
        xbar = IC;
        f_xbar = f(xbar');
        g_xbar = g(xbar');
        [M, N, Gamma, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);
        Pi = Bezier.Pi(c,2,1);
%         
%         A_x_ = [];
%         b_x_ = [];
%         A_u_ = [];
%         b_u_ = [];
%         for m = 1:4
%             I_m = zeros(1,4);
%             I_m(m) = 1;
%             A_tmp = [I_m*H_0'; I_m*H_1'];
%             
%             % input constaints
%             A_u_ = [A_u_; Pi*[eye(2)*A_tmp*D_nT;...
%                 I_m*H_2'*D_nT]];
%             b_u_ = [b_u_; 1+Pi*[xbar; f_xbar]]; % has to match numerator of c.
%             
%             % and state constraints
%             A_x_ = [A_x_; A_x*A_tmp*D_nT];
%             b_x_ = [b_x_; b_x];
%         end
%         A_in = [A_u_; A_x_];
%         b_in = [b_u_; b_x_];
%         
%         A = A_in(:,3:4);
%         b = b_in - A_in(:,1:2)*IC;

        [A_in, b_in] = Bezier.F_G(A_x, b_x, Pi, H, xbar, f_xbar, 2);
        A = A_in; 
        b = b_in;

        
        x0 = IC;
        color = 'g';
                Vert = cddmex('extreme',struct('A',[D(1:2,:); A],'B',[x0;b],'lin',1:2));
        Vert = uniquetol(Vert.V*D(3:4,:)',1e-3,'ByRows',1);
        
%         try
%             Vert = uniquetol(Bezier.Poly.hyp2vert(A,b),1e-3,'ByRows',1);
%         catch E
%             break
%         end
        if ~isempty(Vert)
            if size(Vert,1) > 2
                ind = convhull(Vert);
                
                Vert = Vert(ind,:);
            end
            Vert_tmp = [Vert_tmp; Vert(1:end-1,:)];
        end
        patch(Vert(:,1),Vert(:,2),color,'facealpha',0.1);
    end
    Vert_ = Vert_tmp;
    ind = convhull(Vert_);
    Vert_ = Vert_(ind,:);
        
end
toc
Vert4 = Vert_;

%% Plot results
clf
hold on;
scatter(X0(1),X0(2),50,'filled','k')
patch([b_x(2) b_x(1) -b_x(2) -b_x(1)],[b_x(4) -b_x(3) -b_x(4) b_x(3)],'k','facealpha',0.1)
axis([-b_x(2)-0.1 b_x(1)+0.1 -b_x(4)-1 b_x(3)+1]);
patch(Vert1(:,1),Vert1(:,2),[1 0 0],'facealpha',0.03,'edgecolor',[0.7 0.1 0.1],'linewidth',2);
patch(Vert2(:,1),Vert2(:,2),[0 1 0],'facealpha',0.03,'edgecolor',[0.1 0.7 0.1],'linewidth',2);
patch(Vert3(:,1),Vert3(:,2),[0 0 0],'facealpha',0.03,'edgecolor',[0 0 0],'linewidth',2);
patch(Vert4(:,1),Vert4(:,2),[0 0 1],'facealpha',0.03,'edgecolor',[0.1 0.1 0.7],'linewidth',2);
legend({'','','1 step reachable over N*dt',...
    '1 step B-spline reachable over N*dt, different $$\bar{x}$$',...
    'N step reachable over dt, same $$\bar{x}$$',...
    'N step reachable over dt, different $$\bar{x}$$'},...
    'interpreter','latex','fontsize',15)

%% True reachability
f = @(x) [x(1,:); -sin(x(2,:))];
g = @(x) 1+0*x(1,:);

u_discretization = 50;
X_FL = [];
X_U = [];
tic
x0 = X0;

tic
for i = 1:u_discretization
    for j = 1:u_discretization
%         for k = 1:u_discretization
%             u = @(t) t.^2*(-u_max + (i-1)/(u_discretization-1)*2*u_max)...
%                 +2*t.*(1-t)*(-u_max + (j-1)/(u_discretization-1)*2*u_max)...
%                 +(1-t).^2*(-u_max + (k-1)/(u_discretization-1)*2*u_max);
 u = @(t) t*(-u_max + (i-1)/(u_discretization-1)*2*u_max)...
                +(1-t)*(-u_max + (j-1)/(u_discretization-1)*2*u_max);
            % Input
            [t,x] = ode45(@(t,x) f(x)+g(x)*u(t),[0 horizon_N*dt], x0);
            state_constraint_satisfaction = all(sum(A_x*x'-b_x <= 0,1) == 4);
            input_constraint_satisfaction = 1; % by construction
            if state_constraint_satisfaction && input_constraint_satisfaction
                X_U = [X_U; x(end,:)];
            end
    end
end
toc
tic
for i = 1:u_discretization
    for j = 1:u_discretization
         u = @(t) t*(-u_max + (i-1)/(u_discretization-1)*2*u_max)...
                +(1-t)*(-u_max + (j-1)/(u_discretization-1)*2*u_max);
            % Feedback linearized
            [t,x] = ode45(@(t,x) [0 1; 0 0]*x + [0; 1]*u(t),[0 horizon_N*dt], x0);
            state_constraint_satisfaction = all(sum(A_x*x'-b_x <= 0,1) == 4);
            input_constraint_satisfaction = all(abs(1./g(x').*(-[0 1]*f(x') + u(t)'))<u_max);
            if state_constraint_satisfaction && input_constraint_satisfaction
                X_FL = [X_FL; x(end,:)];
            end
%         end
    end
end
toc
%%
V_FL = Bezier.Poly.conv(X_FL);
V_U= Bezier.Poly.conv(X_U);
patch(V_FL(:,1),V_FL(:,2),[1 1 0],'facealpha',0.03,'edgecolor',[0.7 0.7 0.1],'linewidth',2);
patch(V_U(:,1),V_U(:,2),[0 1 1],'facealpha',0.03,'edgecolor',[0.1 0.7 0.7],'linewidth',2);
l = legend;
l.String{end-1} = 'True Reachable Set of Linear Feedback linearized Inputs';
l.String{end} = 'True Reachable Set of Linear Inputs';
