function [M, N, Gamma, c, M_og] = M_N_Gamma(Lg, Lf, g_xbar, e_bar, K, u_max)

M_og = 1/2*[2*Lg*Lf Lg; Lg 0];
N = [2*Lg*Lf*e_bar + Lf*norm(g_xbar,2) + Lg*norm(K,2)*e_bar; norm(g_xbar,2)+Lg*e_bar];
Gamma = e_bar*(Lg*e_bar + norm(g_xbar,2))*(Lf+norm(K,2));

M = Bezier.Proj_PSD(M_og);

offset = u_max-Gamma;

p1 = [(-N(1) + sqrt(N(1)^2+4*M(1,1)*offset))/(2*M(1,1)) 0];
p2 = [0 (-N(2) + sqrt(N(2)^2+4*M(2,2)*offset))/(2*M(2,2))];
c = [1/p1(1); 1/p2(2)]; % numerator just has to match inequality of c
end