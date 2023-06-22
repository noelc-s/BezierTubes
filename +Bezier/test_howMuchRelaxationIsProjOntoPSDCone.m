%% Description
% A script to see how much of a relaxation the projection onto the PSD cone
% is. To see this, we sample level sets of the quadratic constraint.

% The comparison that we want is arbitrary curves with original bound and
% Bezier curves with new bound, but this is captured by comparing this to
% HJB reachability. Instead, we compare the difference in just the Bezier
% basis, as this is what our method is predicated on.

% Parameters
u_max = 3;
dt = 1;
A_x = [1 0; -1 0; 0 1; 0 -1];
b_x = [0.3; 0.3; 0.6; 0.6];

% Dynamics
f = @(x) -1*sin(x(:,1));
g = @(x) 1+0*x(:,1);
Lf = 1;
Lg = 1; % this is LG_inverse
e_bar = 0;
K = [-1 -1];
% Reference point
x0 = [0; 0];
xbar = [0; 0];
f_xbar = f(xbar');
g_xbar = 1./g(xbar'); % This is g_inverse

[M, N, Gamma, c, M_og] = Bezier.M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max);

% Bezier Matrices
order = 3;
gamma = 2;
H = Bezier.H(order, dt);
D = Bezier.D(gamma,order,dt);
D_nT = inv(D);
Pi = Bezier.Pi(c,2,1);
Z = Bezier.Z(order, dt);

%%
nom = @(sigma) sigma'*M_og*sigma+N'*sigma+Gamma;
ref = @(sigma) sigma'*M*sigma+N'*sigma+Gamma;

density = 100;

[X,Y] = meshgrid(linspace(0,3,density));
Z_nom = zeros(size(X));
Z_ref = zeros(size(X));

Z_k = zeros(size(X));

for i = 1:numel(X)
    Z_nom(i) = nom([X(i) Y(i)]');
    Z_ref(i) = ref([X(i) Y(i)]');
end

clf
subplot(2,1,1)
hold on
[c_nom] = contour(X,Y,Z_nom,[u_max u_max],'color','black');
[c_ref] = contour(X,Y,Z_ref,[u_max u_max],'color','red');

c_nom = c_nom(:,2:end);
c_ref = c_ref(:,2:end);
%%


subplot(2,1,2)
hold on

density = 100;
[X,Y] = meshgrid(linspace(-.3,.3,density),linspace(-.6,.6,density));
for j=1:numel(X)
    
    x1 = [X(j); Y(j)];
    if ~all(A*x1-b <= 0)
        
        xi = D\[x0; x1];
        Xi = [xi H*xi];
        
        % figure(1)
        % scatter(Xi(:,1),Xi(:,2))
        % Z = Bezier.Z(order, dt);
        % tau = linspace(0,1);
        % C = Z(tau)'*Xi;
        % hold on;
        % plot(C(:,1),C(:,2))
        
        q_d_gamma = H^2*xi;
        nom_val = 0;
        ref_val = 0;
        for i = 1:size(Xi,1)
            sigma = [norm(Xi(i,:)'-x0,inf); norm(q_d_gamma(i) - f(x0),inf)];
            nom_val = nom_val + (nom(sigma) - u_max)>0; %should be less than or equal to zero
            ref_val = ref_val + (ref(sigma) - u_max)>0; %should be less than or equal to zero
        end
        subplot(2,1,2)
        if nom_val == 0
            scatter(x1(1),x1(2),30,'filled','k')
        end
        if ref_val == 0
            scatter(x1(1),x1(2),30,'filled','r')
        end
        % subplot(2,1,1) % just plotting the last sigma, that's probably the biggest
        % if nom_val == 0
        %     scatter(sigma(1,end),sigma(2,end),30,'filled','k')
        % end
        % if ref_val == 0
        %     scatter(sigma(1,end),sigma(2,end),30,'filled','r')
        % end
        if mod(j,100)==0
            drawnow
        end
    end
end

% run test_ReachablePend first.
patch(Vert(:,1),Vert(:,2),color,'facealpha',0.1);





