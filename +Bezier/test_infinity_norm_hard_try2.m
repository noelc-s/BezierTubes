clear;

A1 = [1 0; 0 2];
A2 = [1 2];
b1 = [0; 1];
b2 = [-1];
a = [1;0];
c = 10;

Lf = 1;
Lg = 1;
M = [2*Lf*Lg Lg; Lg 0];
M_proj = Bezier.Proj_PSD(M);
N = [1; 2];

offset = c;
p1 = [(-N(1) + sqrt(N(1)^2+4*M_proj(1,1)*offset))/(2*M_proj(1,1)) 0];
p2 = [0 (-N(2) + sqrt(N(2)^2+4*M_proj(2,2)*offset))/(2*M_proj(2,2))];
C = [offset/p1(1); offset/p2(2)]; % numerator just has to match inequality of c

f = @(x) [norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]'*M*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]...
    + N'*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]+a'*x;
f2 = @(x) [norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]'*M_proj*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]...
    + N'*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')];
f3 = @(x) [norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]'*M_proj*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]...
    + N'*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]+a'*x;
f4 = @(x) C'*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]+a'*x;

% try to see when two rows are identical, that's when the derivative is
% undefined
perm = [1 1; 1 -1]; % cols are state space sized.
% perm = [1 1 0; -1 1 0; 0 1 1; 0 -1 1; 1 0 1; 1 0 -1];
x_plot = linspace(-4,4);
for i = 1:size(perm,1)
    line(i,:) = (-(perm(i,:)*A1(:,1))*x_plot + perm(i,:)*b1)/(perm(i,:)*A1(:,2));
end
old_size = size(line,1);
perm = 1;
for i = 1:size(perm,1)
    line(old_size+i,:) = (-(perm(i,:)*A2(:,1))*x_plot + perm(i,:)*b2)/(perm(i,:)*A2(:,2));
end

% l1 = (A1(1,1)*x-b1(1))/A1(1,2);
% l2 = (A1(1,2)*x-b1(2))/A1(2,2);
% l3 = (A2(1,1)*x-b2(1))/A2(1,2);


[X, Y] = meshgrid(linspace(-4,4,100));
Z1 = zeros(size(X));
Z2 = zeros(size(X));
Z3 = zeros(size(X));
Z4 = zeros(size(X));
for i = 1:numel(Z1)
    Z1(i) = f([X(i); Y(i)]);
    Z2(i) = f2([X(i); Y(i)]);
    Z3(i) = f3([X(i); Y(i)]);
    Z4(i) = f4([X(i); Y(i)]);
end

clf;
hold on;
contour(X, Y, Z4,100);
cont1 = contour(X, Y, Z1,[c c]);
cont2 = contour(X, Y, Z2,[c c]);
cont3 = contour(X, Y, Z3,[c c]);

plot(cont1(1,2:end),cont1(2,2:end),'r','linewidth',2)
plot(cont2(1,2:end),cont2(2,2:end),'b','linewidth',2)
plot(cont3(1,2:end),cont3(2,2:end),'k','linewidth',2)
for i = 1:size(line,1)
    plot(x_plot,line(i,:))
end

axis([-4 4 -4 4])

%% something else
% figure(2)
% clf
% hold on
% for k = 1:size(line,1)
%     pts = [x_plot; line(k,:)]
%     for i = 1:size(pts,2)
%         V(i) = f3(pts(:,i));
%     end
%     plot3(pts(1,:),pts(2,:),V)
% end
% plot3(cont2(1,2:end),cont2(2,2:end),10+0*cont2(2,2:end),'b','linewidth',2)
% plot3(cont3(1,2:end),cont3(2,2:end),10+0*cont3(2,2:end),'k','linewidth',2)
% grid on

%
% enumerate each row paris
syms lambda
A = {A1 A1 A2};
b = {b1 b1 b2};
vertices = [];
for aa = 1:size(A,2)
    
%     for i = 1:size(perm,1)
%         line(i,:) = (-(perm(i,:)*A1(:,1))*x_plot + perm(i,:)*b1)/(perm(i,:)*A1(:,2));
%     end
%     old_size = size(line,1);
%     perm = 1;
%     for i = 1:size(perm,1)
%         line(old_size+i,:) = (-(perm(i,:)*A2(:,1))*x_plot + perm(i,:)*b2)/(perm(i,:)*A2(:,2));
%     end
%     perm = [1 0 1; 0 1 1; 1 0 -1; 0 1 -1; -1 0 1; 0 -1 1; -1 0 -1; 0 -1 -1];
%     center1 = pinv(A1)*b1
%     center2 = pinv(A1)*b1
    
    
    z = pinv(A{aa})*b{aa};
    v = [x_plot(1); line(aa,1)];
    x = z+lambda*(v-z);
    lam = [];
    perm = [1 1; 1 -1; -1 1; -1 -1];
    for i = 1:size(A1,1)
        for j = 1:size(A2,1)
            for p = 1:size(perm,1)
                s1_l = perm(p,1)*(A1(i,:)*(v-z));
                s1_c = perm(p,1)*(A1(i,:)*z-b1(i));
                s2_l = perm(p,2)*(A2(j,:)*(v-z));
                s2_c = perm(p,2)*(A2(j,:)*z-b2(j));
                a_ = [s1_l; s2_l]'*M_proj*[s1_l; s2_l];
                b_ = 2*[s1_l; s2_l]'*M_proj*[s1_c; s2_c]+N'*[s1_l; s2_l] + a'*(v-z);
                c_ = [s1_c; s2_c]'*M_proj*[s1_c; s2_c]+N'*[s1_c; s2_c] + a'*z-c;
                l = (-b_+sqrt(b_^2-4*a_*c_))/(2*a_);
                lam = [lam; l];
                l = (-b_-sqrt(b_^2-4*a_*c_))/(2*a_);
                lam = [lam; l];
            end
        end
    end
    % lam = uniquetol(lam,1e-3);
    % [~, ind] = sort(abs(lam));
    % lam = lam(ind(1:2)); % two smallest
    for l = 1:size(lam)
        pt = z+lam(l)*(v-z);
        if abs(f3(pt)-c) < 1e-2
            vertices = [vertices pt];
            scatter(pt(1),pt(2),100,'filled')
        end
    end
end
for i = 1:size(vertices,2)
    level(i) = f4(vertices(:,i));
end
delta = min(level);
cont4 = contour(X, Y, Z4,[delta delta]);
plot(cont4(1,2:end),cont4(2,2:end),'g','linewidth',2)