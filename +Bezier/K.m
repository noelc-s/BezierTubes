function Y = K(m, n)
% stackoverflow "How to compute Commutation matrix in MATLAB"
I = reshape(1:m*n,[m, n]);
I=I';
I=I(:);
Y = eye(m*n);
Y = Y(I,:);
end