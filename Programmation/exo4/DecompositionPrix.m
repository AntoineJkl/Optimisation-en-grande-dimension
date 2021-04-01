function [P_sol,p,k,J,t] = DecompositionPrix(N,P0,a,b,Pmax,parametres,parametres_sousproblemes)
    start = tic;
    %Lien vers les algorithmes 
    addpath('..\Algorithme');

    %Initialisation generale:
    rho = parametres.rho;
    eps = parametres.eps;
    kmax = parametres.kmax;
    
    %Initialisation des sous-problemes:
    option = parametres_sousproblemes.choix;
    rho_sp_uzawa = parametres_sousproblemes.rho_sp_uzawa;
    rho_sp_arrow1 = parametres_sousproblemes.rho_sp_arrow1;
    rho_sp_arrow2 = parametres_sousproblemes.rho_sp_arrow2;
    eps_sp = parametres_sousproblemes.eps_sp;
    kmax_sp = parametres_sousproblemes.kmax_sp;
    mu = zeros(N,1);
    lambda = zeros(N,1);
   
    k = 1;
    P = zeros(N,1);
    P_prec = P + 10;
    P_sol = P;
    p = 0;
    
    while( k <= 2 || ((norm(P - P_prec,2) > eps) && k <= kmax))
        P_prec = P;
        
        %Décomposition:
        for i = 1:N
                %Donnees des sous-problemes:
                A_sp = a(i); b_sp = -p + 2*a(i)*P0(i); C_in = 1; d_in = Pmax(i);
            if option == 1
                %Resolution par Uzawa:
                param_sp_uzawa = struct('rho', rho_sp_uzawa, ...
                        'mu_ini' , mu(i) , ...
                        'lambda_ini' , lambda(i) , ...
                        'eps', eps_sp, ...
                        'kmax', kmax_sp);
                [P(i),lambda(i),mu(i),~] = Uzawa(A_sp,b_sp,0,0,C_in,d_in,param_sp_uzawa);
            else 
                %Resolution par Arrow:
                param_sp_arrow = struct('rho1', rho_sp_arrow1, ...
                        'rho2',rho_sp_arrow2,...
                        'mu_ini' , mu(i) , ...
                        'lambda_ini' , lambda(i) , ...
                        'eps', eps_sp, ...
                        'kmax', kmax_sp);
                [P(i),lambda(i),mu(i),~] = ArrowHurwicz(A_sp,b_sp,0,0,C_in,d_in,param_sp_arrow);
            end
        end
        
        %Coordination:
        p = p + rho*sum(P);
        
        %Recuperation de la solution a l'iteration k:
        P_sol = [P_sol,P];
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
    end
    
    %Calcul de la valeur optimale
    J = sum(a.*(P-P0).^2 + b);
    
    t = toc(start);
end

