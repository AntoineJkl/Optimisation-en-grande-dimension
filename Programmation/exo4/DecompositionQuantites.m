function [P_sol,omega,k,J,t] = DecompositionQuantites(N,P0,a,b,Pmax,parametres,parametres_sousproblemes)    
    start = tic;
    %Pour ajouter les algorithmes 
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
    Mu = zeros(N,1);
    Lambda = zeros(N,1);

    k = 1;    
    P = zeros(N,1);
    P_prec = P + 10;
    %omega = zeros(N,1);
    omega = [1/9*ones(9,1);-1/6*ones(6,1)];
    
    P_sol = P;
       
    while( k <= 2 || ( (norm(P - P_prec) > eps) && k <= kmax))

        P_prec = P;
        
        %Décomposition:
        for i = 1:N
            %Données des sous-problèmes:
            A_sp = a(i); b_sp = 2*a(i)*P0(i); C_in = 1; d_in = Pmax(i); C_eq = 1; d_eq = omega(i);
            if option == 1
                %Resolution par Uzawa:
                param_sp_uzawa = struct('rho', rho_sp_uzawa, ...
                        'mu_ini' , Mu(i) , ...
                        'lambda_ini' , Lambda(i) , ...
                        'eps', eps_sp, ...
                        'kmax', kmax_sp);
                [P(i),Lambda(i),Mu(i),~] = Uzawa(A_sp,b_sp,C_eq,d_eq,C_in,d_in,param_sp_uzawa);
            else 
                %Resolution par Arrow:
                param_sp_arrow = struct('rho1', rho_sp_arrow1, ...
                        'rho2',rho_sp_arrow2,...
                        'mu_ini' , Mu(i) , ...
                        'lambda_ini' , Lambda(i) , ...
                        'eps', eps_sp, ...
                        'kmax', kmax_sp);
                [P(i),Lambda(i),Mu(i),~] = ArrowHurwicz(A_sp,b_sp,C_eq,d_eq,C_in,d_in,param_sp_arrow);
            end    
        end
        
        %Coordination:
        omega = omega + rho*(Lambda - repmat(mean(Lambda),N,1));
        
        P_sol = [P_sol,P];
        
        %Incrementation du nombre d'iterations:
        k = k + 1;

    end
    
    %Calcul de la valeur objectif optimale:
    J = sum(a.*(P-P0).^2 + b);
    
    t = toc(start);
end

