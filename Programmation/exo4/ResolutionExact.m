function [P,J_opt] = ResolutionExact(N,P0,a,b,Pmax)

    A = diag(a);
    b_ = 2*a.*P0;
    c = sum(a.*(P0.^2) + b);

    f = @(x)  x'*A*x - b_'*x + c;

    P_in = zeros(N,1);
    C_in = [];
    d_in = [];

    C_eq = ones(1,N);
    d_eq = 0;
    
    lb = [];
    ub = Pmax;
    options = optimoptions('fmincon','Display','off');
    [P,J_opt] = fmincon(f,P_in,C_in,d_in,C_eq,d_eq,lb,ub,[],options);
end

