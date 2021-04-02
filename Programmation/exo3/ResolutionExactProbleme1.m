function [x_opt,J_opt] = ResolutionExactProbleme1(N,x0,Q,e,Re)

    f = @(x)  .5*x'*Q*x;
    
    C_in = e';
    d_in = Re;

    C_eq = ones(1,N);a
    d_eq = 1;
    
    lb = zeros(1,N);
    ub = [];
    options = optimoptions('fmincon','Display','off');
    [x_opt,J_opt] = fmincon(f,x0,C_in,d_in,C_eq,d_eq,lb,ub,[],options);
end


