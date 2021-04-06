function [P_sol,Multiplicateur_lambda,k,J,t] = DecompositionQuantites(N,P0,a,b,Pmax,parametres,parametres_sousproblemes)    
    start = tic;
    %Pour ajouter les algorithmes : 
    addpath('..\..\Algorithme');
    
    %Initialisation g�n�rale :
    rho = parametres.rho;
    eps = parametres.eps;
    kmax = parametres.kmax;
    
    %Initialisation des sous-probl�mes :
    rho_sp_uzawa = parametres_sousproblemes.rho_sp_uzawa;
    eps_sp = parametres_sousproblemes.eps_sp;
    kmax_sp = parametres_sousproblemes.kmax_sp;
    Mu = zeros(N,1);
    Lambda = zeros(N,1);

    k = 1;    
    P = zeros(N,1);
    P_prec = P + 10;
    omega = zeros(N,1);
    P_sol = [];
    Multiplicateur_lambda = [];
       
    while( k <= 2 || ( (norm(P - P_prec) > eps) && k <= kmax))
        P_prec = P;
        
        disp(['Iterration: ',num2str(k)]);
        
        %D�composition :
        for i = 1:N
            
            %Donn�es des sous-probl�mes :
            A_sp = a(i); b_sp = 2*a(i)*P0(i); C_in = 1; d_in = Pmax(i); C_eq = 1; d_eq = omega(i);
            
            %Resolution par Uzawa :
            param_sp_uzawa = struct('rho', rho_sp_uzawa, ...
                    'mu_ini' , Mu(i) , ...
                    'lambda_ini' , Lambda(i) , ...
                    'eps', eps_sp, ...
                    'kmax', kmax_sp);
            [P(i),Lambda(i),Mu(i),~] = Uzawa(A_sp,b_sp,C_eq,d_eq,C_in,d_in,param_sp_uzawa);
        end
        
        %Coordination :
        omega = omega + rho*(Lambda - repmat(mean(Lambda),N,1));

        %R�cuperation de la solution a chaque it�ration :
        P_sol = [P_sol,P];
        
        %R�cuperation des multiplicateurs � chaque it�ration :
        Multiplicateur_lambda = [Multiplicateur_lambda,mean(Lambda)];
        
        %Incr�mentation du nombre d'it�rations :
        k = k + 1;

    end
    
    %Calcul de la valeur objectif optimale :
    J = sum(a.*(P-P0).^2 + b);
    
    t = toc(start);
end

