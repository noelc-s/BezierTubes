function R = R_n(order_1, order_2)
R_big = eye(order_2+1,order_1+1);
for j = order_1+1:order_2
    R = Bezier.R(j);
    R_big(1:j+1,1:order_1+1) = R*R_big(1:j,1:order_1+1);
end
R = R_big;
end