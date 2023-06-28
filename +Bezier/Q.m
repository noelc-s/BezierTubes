function [Q_store, Q_stack] = Q(segments, order)

M = Bezier.M(order);

Q_stack = zeros(order+1,segments*(order+1));
index = 1;
Q{index} = eye(order+1);
for j = segments-1:-1:0
    
%     if j == 0
%         z = 0;
%     else
        z = j/(j+1);
%     end
    for i = 1:order+1
        Z(i,i) = z^(i-1);
        Z_prime(i,i) = (1-z)^(i-1);
    end
    
    
    Q_1 = M\Z*M; % from zero to z
    Q_2 = fliplr(eye(4))*(M\Z_prime*M)*fliplr(eye(4));
    
    % First segments splitting is iterative
    Q{index+1} = Q_1*Q{index};
    % Second segment gets stored
    Q_store{index} = (Q_2*Q{index})';
    Q_stack(:,(j)*(order+1)+1:(j+1)*(order+1)) = Q_store{index};
    
    index = index+1;
end

% z = 1/segments;
% for i = 1:order+1
%     Z(i,i) = z^(i-1);
% end
% Q_1 = M\Z*M; % from zero to z
% Q_store{index} = Q_1;

% Go from first to last
Q_store = fliplr(Q_store);
Q_stack(:,1:order+1) = Q_store{1};

% Xi_split{j+1} = (Q*Xi')';

end