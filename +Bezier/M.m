function M = M(order)

C = zeros(1,order+1);
for i = 1:order+1
    C(i) = nchoosek(order,i-1);
end
M = diag(C)/abs(pascal(order+1,1)); % forward slash is inverse

%%%%%% alternatively,
% M = pascal(order+1,1);
% for i= 0:order
%     M(i+1,:) = M(i+1,:)*M(end,i+1);
% end
end