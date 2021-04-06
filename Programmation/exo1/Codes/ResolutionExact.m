function [U_opt,J_opt] = ResolutionExact(N)
    [A,b,C] = CreateInstance(N);

    f = @(u)  0.5*u'*A*u - b'*u;

    u0 = zeros(N,1);
    C_in = C;
    d_in = zeros(N,1);

    C_eq = [];
    d_eq = [];
    
    lb = [];
    ub = [];
    
    options = optimoptions('fmincon','Display','off');
    [U_opt,J_opt] = fmincon(f,u0,C_in,d_in,C_eq,d_eq,lb,ub,[],options);

end

