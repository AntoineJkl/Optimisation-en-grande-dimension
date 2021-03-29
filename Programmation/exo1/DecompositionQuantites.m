function [u,omega,k,J,t2] = DecompositionQuantites(N,A,b,C,eps,kmax)
    t1=tic;
    %Pour ajouter les algorithmes (Uzawa et Arrow)
    addpath('..\Algorithme');

    %Initialisation generale:
    k = 1; %Iteration
    u = zeros(N,1); %Solution
    omega = zeros(N,N); %Allocations
    Mu = zeros(N,N); %Multiplicateurs
    
    %Initialisation des sous-problemes:
    rho_sp = 0.1;
    eps_sp = 10^(-4);
    kmax_sp = 10000;
    critere = 0;
    
    %(norm(u - u_prec,2)/norm(u,2) + norm(omega - omega_prec,inf)/norm(omega,inf) > eps)
    
    while( k <= 2 || (~critere && k <= kmax) )
        
        %disp(['iteration: ',num2str(k)]);
        
        u_prec = u;
        omega_prec = omega;
        
        rho=0.1;
        
        %Décomposition:
        for i = 1:N
            A_sp = 1/2*A(i,i);
            b_sp = b(i);
            mu_ini_sp = Mu(:,i);
            C_in_sp = C(:,i);
            d_in = omega(:,i);
            param_sp = struct('rho', rho_sp, ...
                    'mu_ini' , mu_ini_sp , ...
                    'lambda_ini' , 0 , ...
                    'eps', eps_sp, ...
                    'kmax', kmax_sp);
            [u(i),~,Mu(:,i),~] = Uzawa(A_sp,b_sp,0,0,C_in_sp,d_in,param_sp);
        end
        
        %Coordination:
        omega = omega + rho*(Mu - repmat(mean(Mu,2),1,N));
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
        
        %Mise a jour du critere
        %critere = Test_KKT(A,b,C,0,sum(Mu,2),0,0,0,u,10^-1);
        critere = (norm(u - u_prec) < eps);
        %critere = (norm(omega - omega_prec,inf)/norm(omega,inf) < eps);
    end
    
    %Calcul de la valeur optimale
    J = 1/2*u'*A*u - b'*u;
    
    t2=toc(t1);
end

