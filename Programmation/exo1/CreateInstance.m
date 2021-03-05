function [A,b,C] = CreateInstance(N)
    A = eye(N);
    b = (-ones(1,N)).^(2:(N+1));
    b = b';
    C = diag(ones(1,N)) + diag(2*ones(1,N-1),1);
end

