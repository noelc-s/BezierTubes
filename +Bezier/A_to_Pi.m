function [Pi_1, Pi_2, delta] = A_to_Pi(g_xbar, A, b, system, H, x_bar, f_xbar)

Lf = system.L.Lf;
Lg = system.L.Lg;
m = system.m;
n = system.n;
gamma = system.gamma;

P1 = [];
P2 = [];

for i = 1:size(A,1)
    a_ = A(i,:);
    b_ = b(i);
    
    [M, N] = M_N(Lf, Lg, g_xbar, a_(end-1:end));
    
%     a = [a_ 0]';
    A1 = [eye(system.n) zeros(system.n,system.m)];
    A2 = [zeros(system.m,system.n) eye(system.m)];
    b1 = x_bar;
    b2 = f_xbar;
    M_proj = Bezier.Proj_PSD(M);
    
    h = @(x) [norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]'*M_proj*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]...
        + N'*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]+[a_(1:n) zeros(1,m)]*x;
    
    c = MN_to_c(M, N, b_)';
    
    ell = @(x) c*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')] + [a_(1:n) zeros(1,m)]*x;
    
    delta(i) = getDeltaLevelSet(M_proj, N, a_(1:n)', b_, system.n, A1, A2, b1, b2, x_bar, f_xbar, h, ell);
    c = MN_to_c(M, N, delta(i));
    [Pi_1, Pi_2] = c_to_Pi(c, m, gamma);
    P1 = [P1; Pi_1];
    P2 = [P2; Pi_2];
end

Pi_1 = P1;
Pi_2 = P2;
delta=delta';

end

function [M, N] = M_N(Lf, Lg, g_xbar, a)
M = a(2)*1/2*[2*Lg*Lf Lg; Lg 0];
N = [a(2)*Lf*norm(g_xbar,2) + a(1); a(2)*norm(g_xbar,2)];
end

function c = MN_to_c(M, N, delta)
if all((all(M==zeros(2))))
    c = [0; 0];
    return
end
M = Bezier.Proj_PSD(M);

offset = delta;

p1 = [(-N(1) + sqrt(N(1)^2+4*M(1,1)*offset))/(2*M(1,1)) 0];
p2 = [0 (-N(2) + sqrt(N(2)^2+4*M(2,2)*offset))/(2*M(2,2))];
if max(eig(M)) < 1e-5
    c = [0;0];
else
    c = [offset/p1(1); offset/p2(2)]; % numerator just has to match inequality of c
end
end

function [Pi_1, Pi_2] = c_to_Pi(c, m, gamma)
C = [c(1) c(1) -c(1) -c(1);
    c(2) -c(2) c(2) -c(2)];

n = m*gamma;

A = eye(n);
B = eye(m);

Pi_x = kron(A,ones(m,1));
Pi_q = kron(ones(n,1), B);

Pi_1 = [C(1,1)*Pi_x;
    C(1,2)*Pi_x;
    C(1,3)*Pi_x;
    C(1,4)*Pi_x];

Pi_2 = [C(2,1)*Pi_q;
     C(2,2)*Pi_q;
     C(2,3)*Pi_q;
     C(2,4)*Pi_q];
end


function delta = getDeltaLevelSet(M_proj, N,a, b_level_set, n, A1, A2, b1, b2, x_bar, f_xbar, h, ell)

vertices = [];

A = {A1, A2};
b = {b1, b2};

for mat = 1:size(A, 2)
    A_mat = A{mat};
    b_mat = b{mat};
    center = pinv(A_mat)*b_mat;
    if size(A_mat,1) > 1
        inds = nchoosek(1:size(A_mat,1),2);
        inds = [inds; -inds];
    else
        inds = [1 1; 2 2];
    end
    for ind = 1:size(inds,1)
        if inds(ind,1) < 0
            i = -inds(ind,1);
            j = -inds(ind,2);
            sn = -1;
        else
            i = inds(ind,1);
            j = inds(ind,2);
            sn = 1;
        end
        
        if i==j
            % todo: what happens here
            pt = A_mat(1,:)' + 2*center;
        else
            if all((b_mat(i) - sn*b_mat(j)+(A_mat(i,:) - sn*A_mat(j,:))*center)==0)
                pt = pinv(A_mat(i,:) - sn*A_mat(j,:));
            else
                pt = pinv(A_mat(i,:) - sn*A_mat(j,:))*(b_mat(i) - sn*b_mat(j)+(A_mat(i,:) - sn*A_mat(j,:))*center);
            end
        end
        z = center;
        v = pt - center;
        
        lam = [];
        perm = [1 1; 1 -1; -1 1; -1 -1];
        for k = 1:size(A2,1)
            for p = 1:size(perm,1)
                s1_l = perm(p,1)*(A1(i,:)*(v-z)); % i or j would work here
                s1_c = perm(p,1)*(A1(i,:)*z-b1(i));
                s2_l = perm(p,2)*(A2(k,:)*(v-z));
                s2_c = perm(p,2)*(A2(k,:)*z-b2(k));
                a_ = [s1_l; s2_l]'*M_proj*[s1_l; s2_l];
                b_ = 2*[s1_l; s2_l]'*M_proj*[s1_c; s2_c]+N'*[s1_l; s2_l] + a'*(v(1:n)-z(1:n));
                c_ = [s1_c; s2_c]'*M_proj*[s1_c; s2_c]+N'*[s1_c; s2_c] + a'*z(1:n)-b_level_set;
                l = (-b_+sqrt(b_^2-4*a_*c_))/(2*a_);
                lam = [lam; l];
                l = (-b_-sqrt(b_^2-4*a_*c_))/(2*a_);
                lam = [lam; l];
            end
        end
        for l = 1:size(lam)
            pt = z+lam(l)*(v-z);
            if abs(h(pt)-b_level_set) < 1e-3
                vertices = [vertices pt];
%                 scatter(pt(1),pt(2),100,'filled')
            end
        end
    end
end
level = [];
for i = 1:size(vertices,2)
    val = ell(vertices(:,i));
    if val > 0
        level = [level val];
    end
end
delta = min(level);
end
