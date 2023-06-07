function R = R(order, dt)
R = zeros(order+1,order);
for i= 1:order
    R(i,i) = (order+1-i)/order;
    R(i+1,i) = i/order;
end
end