function [P,p,k,J] = DecompositionPrix(N,P0,a,b,Pmax,rho,p_ini,eps,kmax)
    
    %Initialisation generale:
    k = 1;
    P = zeros(N,1);
    p = p_ini;
    
    %Initialisation des sous-problemes:
    rho_sp = 0.2;
    mu_ini_sp = 0;
    lambda_ini_sp = 0;
    eps_sp = 10^(-8);
    kmax_sp = 10000;
    
    while( k <= 2 || ((norm(P - P_prec,2)/norm(P,2) + norm(p - p_prec,2)/norm(p,2) > eps) && k <= kmax))
        disp(['Iteration ',num2str(k)]);
        P_prec = P;
        p_prec = p;
        %Décomposition:
        for i = 1:N
            A = a(i);
            b = -p + 2*a(i)*P0(i);
            C_in = 1;
            d_in = Pmax(i);
            [P(i),~,~,~] = Uzawa(A,b,0,0,C_in,d_in,rho_sp,mu_ini_sp,lambda_ini_sp,eps_sp,kmax_sp);
        end
        %Coordination:
        p = p + rho*sum(P);
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
    end
    
    %Calcul de la valeur optimale
    J = sum(a.*(P-P0).^2 + b);
end

