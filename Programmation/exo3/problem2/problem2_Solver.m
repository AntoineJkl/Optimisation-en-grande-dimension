function [x_opt] = problem2_Solver(e,Q,De, x0)
%PROBLEM2_SOLVER Resolution du probleme 2 par la methode fmincon
    N = length(e);
    if (nargin<4)
        x0 = ones(N, 1)/N; %solution initiale
    end
    
    %objectif
    f = @(x) -e'*x; 
    
    %contraintes
    confuneq = @(x) deal(x'*Q*x - De,...
                         ones(1, N)*x - 1);
    %domaine
    lb = zeros(1,N);
    ub = [];
    
    %resolution par fmincon
    options = optimoptions('fmincon','Display','off','Algorithm','sqp');
    x_opt = fmincon(f,x0,[],[],[],[],lb,ub,confuneq,options);
end


