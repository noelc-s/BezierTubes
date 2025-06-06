function [F, G] = F_G(Ax, bx, H, m, xbar, f_xbar, g_xbar, gamma, Q, Lg, Lf, e_bar, K, u_max)
order = size(H,1);

F = [];
G = [];

k_bar = 0;
e_0 = 0.00;
system.L.L_k = 1;
system.L.L_e = 0.001;
system.L.L_psi = 1;
system.L.L_pi = 1;
system.L.Lf = Lf;
system.L.Lg = Lg;
system.n = 2;
system.m = 1;
system.gamma = 2;
system.order = 3;

C_add = kron(Ax,ones(4*gamma*m*m,1));

for j = 1:size(xbar,2)
    H_vec = Bezier.H_vec(H, m, order-1, gamma, gamma);
%     H_vec2 = Bezier.H_vec(H, m, order-1, gamma, gamma-1);
%     [~,~,~, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar(:,((j-1)*m+1):j*m), e_bar, K, u_max);
%     Pi = Bezier.Pi(c,m, gamma);
    
    % input constraints
    Au = [system.L.L_k*(1+system.L.L_psi) 1+system.L.L_k*system.L.L_e];
    Au = [zeros(size(Au,1), system.n) Au];
    bu = u_max - (norm(k_bar, 2) + e_0);
    [Pi_u_1, Pi_u_2, delta_u] = Bezier.A_to_Pi(g_xbar(:,((j-1)*m+1):j*m), Au, bu, system, H, xbar(:,j), f_xbar(:,j));
    
    % state constraints
    Ax_ = system.L.L_pi*system.L.L_e*sqrt(diag(Ax*Ax'));
    Ax_ = [Ax zeros(size(Ax_,1),1) Ax_];
    bx_ = bx - system.L.L_pi*e_0*sqrt(diag(Ax*Ax'));
    [Pi_x_1, Pi_x_2, delta_x] = Bezier.A_to_Pi(g_xbar(:,((j-1)*m+1):j*m), Ax_, bx_, system, H, xbar(:,j), f_xbar(:,j));
    
    Pi_input_constraints = [Pi_u_1 Pi_u_2];
    Pi_state_constraints = [Pi_x_1 Pi_x_2];
    Pi = [Pi_u_1 Pi_u_2; Pi_x_1+C_add Pi_x_2;];
    h = [kron(delta_u,ones(4*gamma*m*m,1)) + Pi_input_constraints*[xbar(:,j); f_xbar(:,j)];
        kron(delta_x,ones(4*gamma*m*m,1)) + Pi_state_constraints*[xbar(:,j); f_xbar(:,j)]];
    
    F = [F;kron(eye(order),Pi)*kron(Q{j}',kron(eye(m),eye(gamma+1)))*H_vec];
    G = [G;kron(ones(order,1),h)];
    
%     F = [F;kron(eye(order),Pi)*kron(Q{j}',kron(eye(m),eye(gamma+1)))*H_vec; kron(eye(order),Ax)*kron(eye(m),kron(Q{j}',eye(gamma)))*H_vec2];
%     G = [G;kron(ones(order,1),1+Pi*[xbar(:,j); f_xbar(:,j)]); kron(ones(order,1),bx)];
end

end