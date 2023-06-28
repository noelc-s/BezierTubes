function R = R_n(order_1, order_2)
R_big = eye(order_1+1,order_2+1);
for j = order_1+1:order_2
    R = Bezier.R(j);
    R_big(1:order_1+1,1:j+1) = R_big(1:order_1+1,1:j)*R;
end
R = R_big;
end