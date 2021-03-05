function [u,p,k,J] = DecompositionPrediction(N,A,b,C,eps,kmax)
    tic;
    %Pour ajouter les algorithmes (Uzawa et Arrow)
    addpath('..\Algorithme');

    %Initialisation generale:
    k = 1; %Iteration
    u = zeros(N,1); %Solution
    v = zeros(N,1); %Second membre
    p = zeros(N,1); %Prix
    i0 = N; %Sous-probleme choisi
    
    %Parametres de relaxation:
    beta = 0.5;
    gamma = 0.5;
    
    %Initialisation des sous-problemes:
    rho_sp = 0.2;
    eps_sp = 10^(-8);
    kmax_sp = 10000;
    
    while( k <= 2 || ((norm(u - u_prec,2)/norm(u,2) > eps) && k <= kmax))
        u_prec = u;
        
        %Resolution du sous-probleme i0:
        A_sp = 1/2*A(i0,i0);
        b_sp = b(i0);
        mu_ini_sp = zeros(N,1);
        C_in = C(:,i0);
        d_in = v;
        [u(i0),~,Mu,~] = Uzawa(A_sp,b_sp,0,0,C_in,d_in,rho_sp,mu_ini_sp,0,eps_sp,kmax_sp);
        
        %Calcul du prix:
        p = (1-beta)*p + beta*Mu;
        
        %Resolution des autres sous-problemes:
        for i = 1:N
            if i ~= i0
                A_sp = 1/2*A(i,i);
                b_sp = b(i) - C(:,i)'*p;
                [u(i),~,~,~] = Uzawa(A_sp,b_sp,0,0,0,0,rho_sp,0,0,eps_sp,kmax_sp);
            end
        end
        
        %Mise-a-jour de v:
        v = (1-gamma)*v - gamma*C*u;
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
    end
    
    %Calcul de la valeur optimale
    J = 1/2*u'*A*u - b'*u;
    
    toc;
end

