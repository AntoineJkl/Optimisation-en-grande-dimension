function [u,p,k,J,t2] = DecompositionPrediction(N,A,b,C,eps,kmax)
    t1=tic;
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
    Mu = zeros(N,1);
    
    while( k <= 2 || ((norm(u - u_prec,2)/norm(u,2) > eps) && k <= kmax))
        u_prec = u;
        
        %Resolution du sous-probleme i0:
        A_sp = 1/2*A(i0,i0);
        b_sp = b(i0);
        mu_ini_sp = Mu;
        C_in = C(:,i0);
        d_in = v;
        param_0 = struct('rho', rho_sp, ...
                    'mu_ini' , mu_ini_sp , ...
                    'lambda_ini' , 0 , ...
                    'eps', eps_sp, ...
                    'kmax', kmax_sp);
        [u(i0),~,Mu,~] = Uzawa(A_sp,b_sp,0,0,C_in,d_in,param_0);
        
        %Calcul du prix:
        p = (1-beta)*p + beta*Mu;
        
        %Resolution des autres sous-problemes:
        for i = 1:N
            if i ~= i0
                A_sp = 1/2*A(i,i);
                b_sp = b(i) - C(:,i)'*p;
                param_sp = struct('rho', rho_sp, ...
                    'mu_ini' , 0 , ...
                    'lambda_ini' , 0 , ...
                    'eps', eps_sp, ...
                    'kmax', kmax_sp);
                [u(i),~,~,~] = Uzawa(A_sp,b_sp,0,0,0,0,param_sp);
            end
        end
        
        %Mise-a-jour de v:
        v = (1-gamma)*v - gamma*C*u;
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
    end
    
    %Calcul de la valeur optimale
    J = 1/2*u'*A*u - b'*u;
    
    t2=toc(t1);
end

