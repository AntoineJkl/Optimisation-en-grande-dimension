function [u,p,k,J,t2,U_final] = DecompositionPrix(N,A,b,C,rho,eps,kmax,bigU)
    t1=tic;
    %Pour ajouter les algorithmes (Uzawa et Arrow)
    addpath('..\Algorithme');

    %Initialisation generale:
    U_final=[]; %Vecteur des solutions � chaque temps i
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
        %D�composition:
        for i = 1:N
            A_sp = 1/2*A(i,i);
            b_sp = b(i)-C(:,i)'*p;
            param_sp = struct('rho', rho_sp, ...
                    'mu_ini' , 0 , ...
                    'lambda_ini' , 0 , ...
                    'eps', eps_sp, ...
                    'kmax', kmax_sp);
            [u(i),~,~,~] = Uzawa(A_sp,b_sp,0,0,0,0,param_sp);
        end
        %Coordination:
        p = max(0,p + rho*C*u);
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
        
        if bigU
        %Mise � jour des solutions calcul�es
            U_final=[U_final,u];
        end
    end
    
    %Calcul de la valeur optimale
    J = 1/2*u'*A*u - b'*u;
    
    t2=toc(t1);
end

