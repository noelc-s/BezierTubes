function Z = Z(order, dt)
for i = 0:order
    Z{i+1} = @(t) nchoosek(order, i)*(t./dt).^i.*(1-t./dt).^(order-i);
end
Z = @(t) cell2mat(arrayfun(@(idx) Z{idx}(t)', 1:order+1, 'uniform',0));
end