
A1 = [1 0; 0 1];
A2 = [1 2];
b1 = [0; 1];
b2 = [0.1];
a = [0;0];
c = 10;

Lf = 1;
Lg = 1;
M = [2*Lf*Lg Lg; Lg 0];
M_proj = Bezier.Proj_PSD(M);
N = [1; 2];

f = @(x) [norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]'*M*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]...
     + N'*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]+ a'*x;
f2 = @(x) [norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]'*M_proj*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]...
     + N'*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]+ a'*x;
f3 = @(x) [norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]'*M_proj*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')]...
     + N'*[norm(A1*x-b1,'inf'); norm(A2*x-b2,'inf')];

[X, Y] = meshgrid(linspace(-2,2,100));
Z1 = zeros(size(X));
Z2 = zeros(size(X));
Z3 = zeros(size(X));
for i = 1:numel(Z1)
    Z1(i) = f([X(i); Y(i)]);
    Z2(i) = f2([X(i); Y(i)]);
    Z3(i) = f3([X(i); Y(i)]);
end

% surf(X, Y, Z)

%%% from delta to level set
delta = c;
p1 = [(-N(1) + sqrt(N(1)^2+4*M_proj(1,1)*delta))/(2*M_proj(1,1)) 0];
p2 = [0 (-N(2) + sqrt(N(2)^2+4*M_proj(2,2)*delta))/(2*M_proj(2,2))];
if max(eig(M_proj)) < 1e-5
    const = [0;0];
else
    const = [delta/p1(1); delta/p2(2)]; % numerator just has to match inequality of c
end
C = [const(1) const(1) -const(1) -const(1);
    const(2) -const(2) const(2) -const(2)];
n = 2;
m=1;
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
Pi = [Pi_1*A1+Pi_2*A2];
h = delta+Pi_1*b1+Pi_2*b2;
Vert = cddmex('extreme',struct('A',Pi,'B',h));
Vert = Bezier.Poly.conv(Vert.V);
if ~isempty(Vert)
    ind = convhull(Vert);
    Vert = Vert(ind,:);
end
%%%

perm = [1 1 1; 1 1 -1; 1 -1 1; -1 1 1; 1 -1 -1; -1 1 -1; -1 -1 1; -1 -1 -1];
center = pinv([A1; A2])*[b1; b2];
v = Vert';

for i = 1:size(v,2)
        v_i = v(:,i);
        a_ = [norm(A1*v_i,'inf'); norm(A2*v_i,'inf')]'*M_proj*[norm(A1*v_i,'inf'); norm(A2*v_i,'inf')];
        b_ = a'*v_i + N'*[norm(A1*v_i,'inf'); norm(A2*v_i,'inf')];
        c_ = -c;
        l_star = (-b_+sqrt(b_^2-4*a_*c_))/(2*a_);
        lambda(i) = l_star;
%         lambda(i) = 2;
end

%%% from delta to level set
for i = 1:size(v,2)
    val(i) = f3(center+lambda(i)*v(:,i));
end
delta = min(val);
p1 = [(-N(1) + sqrt(N(1)^2+4*M_proj(1,1)*delta))/(2*M_proj(1,1)) 0];
p2 = [0 (-N(2) + sqrt(N(2)^2+4*M_proj(2,2)*delta))/(2*M_proj(2,2))];
if max(eig(M_proj)) < 1e-5
    const = [0;0];
else
    const = [delta/p1(1); delta/p2(2)]; % numerator just has to match inequality of c
end
C = [const(1) const(1) -const(1) -const(1);
    const(2) -const(2) const(2) -const(2)];
n = 2;
m=1;
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
Pi = [Pi_1*A1+Pi_2*A2];
h = delta+Pi_1*b1+Pi_2*b2;
Vert = cddmex('extreme',struct('A',Pi,'B',h));
Vert = Bezier.Poly.conv(Vert.V);
if ~isempty(Vert)
    ind = convhull(Vert);
    Vert = Vert(ind,:);
end
%%%

cont1 = contour(X, Y, Z1,[c c]);
cont2 = contour(X, Y, Z2,[c c]);
cont3 = contour(X, Y, Z3,.7*[delta delta]);
clf
hold on
plot(cont1(1,2:end),cont1(2,2:end),'r')
plot(cont2(1,2:end),cont2(2,2:end),'b')
plot(cont3(1,2:end),cont3(2,2:end),'k')
% patch(Vert(:,1),Vert(:,2),[1 0 0],'facealpha',0.03,'edgecolor',[0.7 0.1 0.1],'linewidth',2);
for i = 1:size(v,2)
%     line([center(1) center(1)+lambda(i)*v(1,i)],...
%         [center(2) center(2)+lambda(i)*v(2,i)]);
    line([0 lambda(i)*v(1,i)],...
        [0 lambda(i)*v(2,i)]);
end
axis([-2 2 -2 2])

% delta = min(vecnorm(lambda'.*v,'inf',2))^2;
% 
% N = zeros(size(X));
% for i = 1:numel(Z)
%     N(i) = norm([X(i); Y(i)],'inf')^2;
% end
% cont_lb = contour(X, Y, N,[delta delta]);
% 
% clf
% hold on
% plot(cont(1,2:end),cont(2,2:end),'r')
% plot(cont_lb(1,2:end),cont_lb(2,2:end),'b')
% axis equal