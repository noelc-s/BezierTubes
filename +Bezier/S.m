
function S = S(order, dt)
S = zeros(order+1,order);
for i= 1:order
    S(i,i) = -order;
    S(i+1,i) = order;
end
end

