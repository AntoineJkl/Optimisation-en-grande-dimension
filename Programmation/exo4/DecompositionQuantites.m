function [P_sol,omega,k,J] = DecompositionQuantites(N,P0,a,b,Pmax,eps,kmax)    
    tic;
    %Pour ajouter les algorithmes 
    addpath('..\Algorithme');

    %Initialisation generale:
    k = 1;
    rho = 0.2;
    P = zeros(N,1);
    Lambda = zeros(N,1);
    Mu = zeros(N,1);
    omega = zeros(N,1);
    
    P_sol = P;
    
    %Initialisation des sous-problemes:
    rho_sp = 0.005;
    eps_sp = 10^(-3);
    kmax_sp = 5000;
    
    while( k <= 2 || ( (norm(P - P_prec) > eps) && k <= kmax))
        
        if mod(k,1) ==  0
            disp(['Iteration: ',num2str(k)]);
        end
        
        %Precedent:
        P_prec = P;
        
        %Décomposition:
        for i = 1:N
            A_sp = a(i);
            b_sp = 2*a(i)*P0(i);
            C_in = 1;
            d_in = Pmax(i);
            C_eq = 1;
            d_eq = omega(i);
            param_sp = struct('rho', rho_sp, ...
                    'mu_ini' , Mu(i) , ...
                    'lambda_ini' , Lambda(i) , ...
                    'eps', eps_sp, ...
                    'kmax', kmax_sp);
            [P(i),Lambda(i),Mu(i),~] = Uzawa(A_sp,b_sp,C_eq,d_eq,C_in,d_in,param_sp);
            
            %[P(i),Lambda(i),~,~] = ArrowHurwicz(A_sp,b_sp,C_eq,d_eq,C_in,d_in,rho_sp,rho_sp,mu_ini_sp,lambda_ini_sp,eps_sp,kmax_sp);
            %P(i) = omega(i); Lambda(i) = -2*a(i)*P(i) + 2*a(i)*P0(i);
        end
        
        %Coordination:
        omega = omega + rho*(Lambda - repmat(mean(Lambda),N,1));
        
        P_sol = [P_sol,P];
        
        %Incrementation du nombre d'iterations:
        k = k + 1;

    end
    
    %Calcul de la valeur objectif optimale:
    J = sum(a.*(P-P0).^2 + b);
    
    toc;
end

