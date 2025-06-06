
a = [.3; .2];
c = .1;

f = @(x) norm(x,'inf')^2 + a'*x;

[X, Y] = meshgrid(linspace(-2,2,200));
Z = zeros(size(X));
for i = 1:numel(Z)
    Z(i) = f([X(i); Y(i)]);
end

cont = contour(X, Y, Z,[c c]);


v = [1 1; 1 -1; -1 1; -1 -1];
for i = 1:size(v,1)
    v_i = v(i,:)';
    a_ = norm(v_i,'inf')^2;
    b_ = a'*v_i;
    c_ = -c;
    l_star = (-b_+sqrt(b_^2-4*a_*c_))/(2*a_);
    lambda(i) = l_star;
%     lambda(i) = c / (norm(v_i,'inf')^2 + a'*v_i);
end
delta = min(vecnorm(lambda'.*v,'inf',2))^2;

N = zeros(size(X));
for i = 1:numel(Z)
    N(i) = norm([X(i); Y(i)],'inf')^2;
end
cont_lb = contour(X, Y, N,[delta delta]);

clf
hold on
plot(cont(1,2:end),cont(2,2:end),'r')
plot(cont_lb(1,2:end),cont_lb(2,2:end),'b')
axis equal