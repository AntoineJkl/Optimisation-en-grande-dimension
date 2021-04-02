function [x_opt,J_opt] = ResolutionExactProbleme2(N,x0,e,Q,De)

    f = @(x) -e'*x;
    
    confuneq = @(x) deal( x'*Q*x - De, sum(x) - 1 );
    
    C_in = [];
    d_in = [];

    C_eq = [];
    d_eq = [];
    
    lb = zeros(1,N);
    ub = [];
    options = optimoptions('fmincon','Display','off','Algorithm','sqp');
    [x_opt,J_opt] = fmincon(f,x0,C_in,d_in,C_eq,d_eq,lb,ub,confuneq,options);
end


