




function S = S(order, dt)
S = zeros(order,order+1);
for i= 1:order
    S(i,i) = -order;
    S(i,i+1) = order;
end
end

