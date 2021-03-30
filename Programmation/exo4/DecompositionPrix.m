function [P_sol,p,k,J] = DecompositionPrix(N,P0,a,b,Pmax,rho,eps,kmax)
    tic;
    %Pour ajouter les algorithmes 
    addpath('..\Algorithme');

    %Initialisation generale:
    k = 1;
    P = zeros(N,1);
    p = 0;
    
    P_sol = P;
    
    %Initialisation des sous-problemes:
    rho_sp = 0.3;
    mu = zeros(N,1);
    lambda = zeros(N,1);
    eps_sp = 10^(-5);
    kmax_sp = 10000;
    
    while( k <= 2 || ((norm(P - P_prec,2) > eps) && k <= kmax))
        P_prec = P;
        %Décomposition:
        for i = 1:N
            A_sp = a(i);
            b_sp = -p + 2*a(i)*P0(i);
            C_in = 1;
            d_in = Pmax(i);
            param_sp = struct('rho', rho_sp, ...
                    'mu_ini' , mu(i) , ...
                    'lambda_ini' , lambda(i) , ...
                    'eps', eps_sp, ...
                    'kmax', kmax_sp);
            [P(i),lambda(i),mu(i),~] = Uzawa(A_sp,b_sp,0,0,C_in,d_in,param_sp);
        end
        %Coordination:
        p = p + rho*sum(P);
        
        P_sol = [P_sol,P];
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
    end
    
    %Calcul de la valeur optimale
    J = sum(a.*(P-P0).^2 + b);
    
    toc;
end

