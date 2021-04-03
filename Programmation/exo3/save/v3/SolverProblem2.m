function [x_opt] = SolverProblem2(e,Q,De, x0)

    N = length(e);
    if (nargin<4)
        x0 = ones(N, 1)/N;
    end
    
    f = @(x) -e'*x;
    
    confuneq = @(x) deal( x'*Q*x - De, ones(1, N)*x - 1 );
    
    C_in = [];
    d_in = [];

    C_eq = [];
    d_eq = [];
    
    lb = zeros(1,N);
    ub = [];
    options = optimoptions('fmincon','Display','off','Algorithm','sqp');
    x_opt = fmincon(f,x0,C_in,d_in,C_eq,d_eq,lb,ub,confuneq,options);
end


