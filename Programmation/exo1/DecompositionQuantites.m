function [u,omega,k,J] = DecompositionQuantites(N,A,b,C,eps,kmax)
    tic;
    %Pour ajouter les algorithmes (Uzawa et Arrow)
    addpath('..\Algorithme');

    %Initialisation generale:
    k = 1; %Iteration
    u = zeros(N,1); %Solution
    omega = zeros(N,N); %Allocations
    Mu = zeros(N,N); %Multiplicateurs
    
    %Initialisation des sous-problemes:
    rho_sp = 0.2;
    eps_sp = 10^(-8);
    kmax_sp = 10000;
    
    while( k <= 2 || ((norm(u - u_prec,2)/norm(u,2) + norm(omega - omega_prec,inf)/norm(omega,inf) > eps) && k <= kmax))
        u_prec = u;
        omega_prec = omega;
        
        %Pas pour la coordination:
        rho = 1/k^2;
        
        %Décomposition:
        for i = 1:N
            A_sp = 1/2*A(i,i);
            b_sp = b(i);
            mu_ini_sp = zeros(N,1);
            C_in = C(:,i);
            d_in = omega(:,i);
            [u(i),~,Mu(:,i),~] = Uzawa(A_sp,b_sp,0,0,C_in,d_in,rho_sp,mu_ini_sp,0,eps_sp,kmax_sp);
        end
        
        %Coordination:
        omega = omega + rho*(Mu - repmat(mean(Mu,2),1,N));
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
    end
    
    %Calcul de la valeur optimale
    J = 1/2*u'*A*u - b'*u;
    
    toc;
end

