function [x_opt,J_opt] = ResolutionExactProbleme1(Q,e,Re,x0)
    
    N = length(e);
    if (nargin<4)
        x0 = ones(N, 1)/N;
    end
    
    f = @(x)  .5*x'*Q*x;
    
    C_in = -e';
    d_in = -Re;

    C_eq = ones(1,N);
    d_eq = 1;
    
    lb = zeros(1,N);
    ub = [];
    options = optimoptions('fmincon','Display','off');
    [x_opt,J_opt] = fmincon(f,x0,C_in,d_in,C_eq,d_eq,lb,ub,[],options);
end


