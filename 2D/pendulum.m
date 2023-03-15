init
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

% FL
syms v
u = LgLfy\(-Lf2y + v);
FL_u = matlabFunction(u,'Vars',[x1, x2, v]);

% CLF outputs
eta = matlabFunction([y; Lfy],'Vars',x_sym);
clear x1 x2

%%

dt = 1;
x0 = [0 0];
u_max = 1;

Lf = 1;
Lg = 1;

x_bar = [0;0];
f_xbar = f_func(x_bar(1), x_bar(2));
g_xbar = g_func(x_bar(1), x_bar(2));

A_x = [0 1; 0 -1; 1 0; -1 0];
b_x = [1; 1; 1; 1];

M_og = 1/2*[2*Lg*Lf Lg; Lg 0];
N = [Lf*norm(g_xbar,2); norm(g_xbar,2)];
M = projectOntoSemidefiniteCone(M_og);

f = figure(1);
% f.Position = [0 0 650 650];

order = 3;
[H, D_nT] = Poly.getBezMatrices(order, dt);

% Inputs
%%% Constraint on first point:
%%% Wait a second, this IS the set.
H_0 = H^0;
H_1 = H^1;
H_2 = H^2;

clf;
hold on;
[X,Y] = meshgrid(linspace(-10,10));
ind_u = ones(size(X));
A_u = [];
b_u = [];

p1 = [(-N(1) + sqrt(N(1)^2+4*M(1,1)*u_max))/(2*M(1,1)) 0];
p2 = [0 (-N(2) + sqrt(N(2)^2+4*M(2,2)*u_max))/(2*M(2,2))];
c_ = [1/p1(1); 1/p2(2)];


% x0 such that x0,x1,x2,x3 satisfy the input constraints
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

s_u = size(A_u,1);
s_x = size(A_x,1);

A_in = [A_u; A_x_];
b_in = [b_u; b_x_];

Vert = Poly.conv(Poly.hyp2vert(A_in,b_in));
patch(Vert(:,1),Vert(:,2),'b','facealpha',0.1);

%%% Linear system
[A_lin, b_lin] = Poly.forwardReachable(H, Df_func(x_bar(1),x_bar(2)),...
    g_func(x_bar(1), x_bar(2)), A_x, b_x,u_max,D_nT,x0);
Vert_lin = Poly.conv(Poly.hyp2vert(A_lin,b_lin));
patch(Vert_lin(:,1),Vert_lin(:,2),'k','facealpha',0.1);

%%% Plot interpolations of the edges
% x1 = x0;
tau = linspace(0,1);
for i = 1:size(Vert,1)-1
    for lambda = 0:.25:1
        X = (1-lambda)*Vert(i,1) + lambda*Vert(i+1,1);
        Y = (1-lambda)*Vert(i,2) + lambda*Vert(i+1,2);
        
        x1 = [X Y];
        a = x0(1);
        b = 1/3*(dt*x0(2)+3*x0(1));
        c = 1/3*(-dt*x1(2)+3*x1(1));
        d = x1(1);
        pBez = @(t) a*(1-t).^3 + 3*b*t.*(1-t).^2 + 3*c*t.^2.*(1-t) + d*t.^3;
        vBez = @(t) (3*d*t.^2 - 3*c*t.^2 - 3*a*(t - 1).^2 + 3*b*(t - 1).^2 - 6*c*t.*(t - 1) + 3*b*t.*(2*t - 2))*1/dt;
        
        scatter([a b c d], [a b c d]*H,30,'b','filled')
        p3=plot(pBez(tau),vBez(tau),'color','b');
        p4=plot([a b c d],[a b c d]*H,'color','b','linestyle','--');
    end
end
axis([-1 1 -1 1])

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