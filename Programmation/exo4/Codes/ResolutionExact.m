function [P,J_opt,lambda] = ResolutionExact(N,P0,a,b,Pmax)

    %Fonction objectif :
    A = diag(a);
    b_ = 2*a.*P0;
    c = sum(a.*(P0.^2) + b);
    f = @(x)  x'*A*x - b_'*x + c;
    
    %Valeur initiale :
    P_in = zeros(N,1);
    
    %Contraintes d'inégalité :
    C_in = [];
    d_in = [];

    %Contraintes d'égalité :
    C_eq = ones(1,N);
    d_eq = 0;
    
    %Bornes sur la solution :
    lb = [];
    ub = Pmax;
    
    %Résolution :
    options = optimoptions('fmincon','Display','off');
    [P,J_opt,~,~,multiplicateur] = fmincon(f,P_in,C_in,d_in,C_eq,d_eq,lb,ub,[],options);
    lambda = multiplicateur.eqlin;
end

