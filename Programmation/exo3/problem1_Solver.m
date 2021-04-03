function [x_opt] = problem1_Solver(Q,e,Re,x0)
%PROBLEM1_SOLVER Resolution du probleme 1 par la methode fmincon
    
    N = length(e);
    if (nargin<4)
        x0 = ones(N, 1)/N;
    end
    
    %objectif
    f = @(x)  .5*x'*Q*x;
    
    %contraintes d'inegalite
    C_in = -e';
    d_in = -Re;
    
    %contraintes d'egalite
    C_eq = ones(1,N);
    d_eq = 1;
    
    %bornes du domaine
    lb = zeros(1,N);
    ub = [];
    
    %resolution par fmincon
    options = optimoptions('fmincon','Display','off');
    x_opt = fmincon(f,x0,C_in,d_in,C_eq,d_eq,lb,ub,[],options);
end


