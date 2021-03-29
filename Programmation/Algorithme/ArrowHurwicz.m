function [U,Lambda,Mu,k] = ArrowHurwicz(A,b,C_eq,d_eq,C_in,d_in,rho1,rho2,mu_ini,lambda_ini,eps,kmax,U_ini)
    
    %Arguments d'entree
    switch nargin 
        case 8
            mu_ini = zeros(size(C_in,1),1);
            lambda_ini = zeros(size(C_eq,1),1);            
            kmax = 10000;
            eps = 10^(-5);
        case 10
            kmax = 10000;
            eps = 10^(-5);
        case 11
            kmax = 10000;
    end
    
    %Initialisation
    k = 1;
    U = U_ini;
    Mu = mu_ini;
    Lambda = lambda_ini;
    
    %Resolution probleme
    while(k <= 2 || (norm( U - prec,2)/norm(U,2) > eps) && (k <= kmax) )
        %Stockage de U precedent
        prec = U;
        %Récupération de U (Formule de projection)
        U = U - rho1*(2*A*U - b + C_eq'*Lambda + C_in'*Mu);
        %Mise a jour des multiplicateurs de Lagrange
        Lambda = Lambda + rho2*(C_eq*U - d_eq);
        Mu = max(0,Mu + rho2*(C_in*U - d_in));
        %Incrementation du nombre d'iterations
        k = k + 1;
    end
end

