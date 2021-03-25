function [P,omega,k,J] = DecompositionQuantites(N,P0,a,b,Pmax,eps,kmax)    
    tic;
    %Pour ajouter les algorithmes 
    addpath('..\Algorithme');

    %Initialisation generale:
    k = 1;
    P = zeros(N,1);
    Lambda = zeros(N,1);
    %omega = zeros(N,1);
    %omega = [15/18*ones(1,9),-15/12*ones(1,6)]';
    omega = [15/8;-15/22*ones(5,1);15/8*ones(3,1);-15/22*ones(6,1)];
    
    %Initialisation des sous-problemes:
    rho_sp = 0.1;
    eps_sp = 10^(-5);
    kmax_sp = 50000;
    
    %Critere d'arret:
    critere = 1;
    
    while( k <= 2 || (critere && k <= kmax))
        if(abs(sum(P))> 10000)
            break;
        end
        disp(['Iteration: ',num2str(k)]);
        
        %Mise a jour du pas:
        rho = 0.005;
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
            [P(i),Lambda(i),~,~] = Uzawa(A_sp,b_sp,C_eq,d_eq,C_in,d_in,rho_sp,mu_ini_sp,lambda_ini_sp,eps_sp,kmax_sp);
            %[P(i),Lambda(i),~,~] = ArrowHurwicz(A_sp,b_sp,C_eq,d_eq,C_in,d_in,rho_sp,rho_sp,mu_ini_sp,lambda_ini_sp,eps_sp,kmax_sp);
            %P(i) = omega(i); Lambda(i) = -2*a(i)*P(i) + 2*a(i)*P0(i);
        end
        %Coordination:
        omega = omega + rho*(Lambda - repmat(1/N*sum(Lambda),N,1));
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
        
        %Mise a jour du critere:
        critere = (norm(P - P_prec,2)/norm(P + 10^(-6),2) + norm(omega - omega_prec,2)/norm(omega + 10^(-6),2) > eps);
    end
    
    %Calcul de la valeur optimale
    J = sum(a.*(P-P0).^2 + b);
    
    toc;
end

