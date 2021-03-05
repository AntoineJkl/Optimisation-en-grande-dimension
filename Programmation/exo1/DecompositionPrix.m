function [u,p,k,J] = DecompositionPrix(N,A,b,C,rho,eps,kmax)
    tic;
    %Pour ajouter les algorithmes (Uzawa et Arrow)
    addpath('..\Algorithme');

    %Initialisation generale:
    k = 1; %Iteration
    u = zeros(N,1); %Solution
    p = zeros(N,1); %Prix
    
    %Initialisation des sous-problemes:
    rho_sp = 0.1;
    eps_sp = 10^(-5);
    kmax_sp = 10000;
    
    while( k <= 2 || ((norm(u - u_prec,2)/norm(u,2) + norm(p - p_prec,2)/norm(p,2) > eps) && k <= kmax))
        u_prec = u;
        p_prec = p;
        %Décomposition:
        for i = 1:N
            A_sp = 1/2*A(i,i);
            b_sp = b(i)-C(:,i)'*p;
            [u(i),~,~,~] = Uzawa(A_sp,b_sp,0,0,0,0,rho_sp,0,0,eps_sp,kmax_sp);
        end
        %Coordination:
        p = max(0,p + rho*C*u);
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
    end
    
    %Calcul de la valeur optimale
    J = 1/2*u'*A*u - b'*u;
    
    toc;
end

