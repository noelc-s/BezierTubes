function [F, G] = F_G(Ax, bx, H, m, xbar, f_xbar, g_xbar, gamma, varargin)

order = size(H,1);

use_Q = false;

if ~isempty(varargin)
    Q = varargin{1};
    Lg = varargin{2};
    Lf = varargin{3};
    e_bar = varargin{4};
    K = varargin{5};
    u_max = varargin{6};
%     use_Q = true;
else
    Q = {eye(order)};
end

F = [];
G = [];

for j = 1:size(xbar,2)
    H_vec = Bezier.H_vec(H, m, order-1, gamma, gamma);
    H_vec2 = Bezier.H_vec(H, m, order-1, gamma, gamma-1);
    [~,~,~, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar(:,((j-1)*m+1):j*m), e_bar, K, u_max);
    Pi = Bezier.Pi(c,m, gamma);
    F = [F;kron(eye(order),Pi)*kron(Q{j}',kron(eye(m),eye(gamma+1)))*H_vec; kron(eye(order),Ax)*kron(eye(m),kron(Q{j}',eye(gamma)))*H_vec2];
    G = [G;kron(ones(order,1),1+Pi*[xbar(:,j); f_xbar(:,j)]); kron(ones(order,1),bx)];
end

end