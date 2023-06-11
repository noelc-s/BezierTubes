function M_ = Proj_PSD(M)
% project onto semidefinite cone
[evec,eval] = eig(M);
M_ = eps*ones(size(M,1));
for i = 1:size(M,1)
    if eval(i,i) > 0
        M_ = M_ + evec(:,i)*evec(:,i)'*eval(i,i);
    end
end
end