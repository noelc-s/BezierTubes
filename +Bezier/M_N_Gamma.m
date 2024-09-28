function [M, N, Gamma, c, M_og] = M_N_Gamma(L, g_xbar, e_bar, u_max, a, b)

% M_og = 1/2*[2*Lg*Lf Lg; Lg 0];
% N = [2*Lg*Lf*e_bar + Lf*norm(g_xbar,2) + Lg*norm(K,2)*e_bar; norm(g_xbar,2)+Lg*e_bar];
% Gamma = e_bar*(Lg*e_bar + norm(g_xbar,2))*(Lf+norm(K,2));

k_bar = 0;
L_k = L.L_k;
L_e = L.L_e;
L_psi = L.L_psi;
L_pi = L.L_pi;
Lf = L.L_f;
Lg = L.L_g;

% for input constraints
% a = [L_k*(1+L_psi) L_k*L_e];
% b = u_max - (norm(k_bar, 2) + e_bar)

% for state constraints
% A = L_pi*L_e*sqrt(diag(C*C'))
% A = [zeros(size(A,1),1) A]
% b = d - L_pi*e_0*sqrt(diag(C*C'))

M_og = a(2)*1/2*[2*Lg*Lf Lg; Lg 0];
N = [a(2)*Lf*norm(g_xbar,2) - a(1); a(2)*norm(g_xbar,2)];

M = Bezier.Proj_PSD(M_og);

offset = b;

p1 = [(-N(1) + sqrt(N(1)^2+4*M(1,1)*offset))/(2*M(1,1)) 0];
p2 = [0 (-N(2) + sqrt(N(2)^2+4*M(2,2)*offset))/(2*M(2,2))];
if max(eig(M)) < 1e-8
    c = [0;1./offset];
else
    c = [1/p1(1); 1/p2(2)]; % numerator just has to match inequality of c
end
end