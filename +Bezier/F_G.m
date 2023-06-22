function [F, G] = F_G(Ax, bx, H, xbar, f_xbar, g_xbar, gamma, varargin)

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

% input constraints

A_u = [];
b_u = [];
A_x = [];
b_x = [];

for j = 1:size(xbar,2)
for i = 1:order
    I_m = zeros(1,order);
    I_m(i) = 1;
    I_m = I_m*Q{j};
%     if use_Q
        [~,~,~, c] = Bezier.M_N_Gamma(Lg, Lf, g_xbar(:,j), e_bar, K, u_max);
        Pi = Bezier.Pi(c,2,1); % n,m
%     end
    A_u = [A_u; Pi*kron(eye(gamma+1),I_m)];
    b_u = [b_u; 1+Pi*[xbar(:,j); f_xbar(:,j)]];
    
    A_x = [A_x; Ax*kron(eye(gamma),I_m)];
    b_x = [b_x; bx];
end
end
H_ = [];
% Q_i times H
for i = 1:gamma
    H_ = [H_; H^(i-1)];
end

F = [A_u*[H_; H^gamma]; A_x*H_];
G = [b_u; b_x];

% A_u = [ Pi*kron(eye(gamma+1),[0 1 0 0]); Pi*kron(eye(3),[0 0 1 0]); Pi*kron(eye(3),[0 0 0 1])]*;
% b_u = kron(ones(order,1),1+Pi*[xbar; f_xbar]);    

% state constraints
% A_x = [A_x*kron(eye(2),[1 0 0 0]); A_x*kron(eye(2),[0 1 0 0]); A_x*kron(eye(2),[0 0 1 0]); A_x*kron(eye(2),[0 0 0 1])]*[H^0; H^1];
% b_x = kron(ones(4,1),b_x);

end