function [P,omega,k,J] = DecompositionQuantites(N,P0,a,b,Pmax,eps,kmax)
    tic;
    %Pour ajouter les algorithmes 
    addpath('..\Algorithme');

    %Initialisation generale:
    k = 1;
    P = zeros(N,1);
    Lambda = zeros(N,1);
    omega = zeros(N,1);
    
    %Initialisation des sous-problemes:
    rho_sp = 0.05;
    eps_sp = 10^(-3);
    kmax_sp = 50000;
    
    while( k <= 2 || ((1 + norm(P - P_prec,2)/norm(P + 10^(-6),2) + norm(omega - omega_prec,2)/norm(omega + 10^(-6),2) > eps) && k <= kmax))
        
        disp(['Iteration: ',num2str(k)]);
        
        %Mise a jour du pas:
        rho = 1/k^2;
        P_prec = P;
        omega_prec = omega;
        
        %Décomposition:
        for i = 1:N
            A_sp = a(i);
            b_sp = 2*a(i)*P0(i);
            C_in = 1;
            d_in = Pmax(i);
            C_eq = 1;
            d_eq = omega(i);
            mu_ini_sp = 0;
            lambda_ini_sp = 0;
            [P(i),Lambda(i),~,n] = Uzawa(A_sp,b_sp,C_eq,d_eq,C_in,d_in,rho_sp,mu_ini_sp,lambda_ini_sp,eps_sp,kmax_sp);
        end
        %Coordination:
        omega = omega + rho*(Lambda - repmat(1/N*sum(Lambda),N,1));
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
    end
    
    %Calcul de la valeur optimale
    J = sum(a.*(P-P0).^2 + b);
    
    toc;
end

