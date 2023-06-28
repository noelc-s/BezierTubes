function R = R(order, dt)
R = zeros(order,order+1);
for i= 1:order
    R(i,i) = (order+1-i)/order;
    R(i,i+1) = i/order;
end
end