function [P_sol,Prix,k,J,t] = DecompositionPrix(N,P0,a,b,Pmax,parametres,parametres_sousproblemes)
    start = tic;
    %Pour ajouter les algorithmes :
    addpath('..\..\Algorithme');

    %Initialisation generale :
    rho = parametres.rho;
    eps = parametres.eps;
    kmax = parametres.kmax;
    
    %Initialisation des sous-probl�mes :
    rho_sp_uzawa = parametres_sousproblemes.rho_sp_uzawa;
    eps_sp = parametres_sousproblemes.eps_sp;
    kmax_sp = parametres_sousproblemes.kmax_sp;
    mu = zeros(N,1);
   
    k = 1;
    P = zeros(N,1);
    P_prec = P + 10;
    p = 0;
    P_sol = P;
    Prix = p;    
    
    while( k <= 2 || ((norm(P - P_prec,2) > eps) && k <= kmax))
        P_prec = P;
        
        disp(['Iterration: ',num2str(k)]);
        
        %D�composition :
        for i = 1:N
            
                %Donn�es des sous-probl�mes :
                A_sp = a(i); b_sp = -p + 2*a(i)*P0(i); C_in = 1; d_in = Pmax(i);
                
                %R�solution par Uzawa :
                param_sp_uzawa = struct('rho', rho_sp_uzawa, ...
                        'mu_ini' , mu(i) , ...
                        'eps', eps_sp, ...
                        'kmax', kmax_sp);
                [P(i),~,mu(i),~] = Uzawa(A_sp,b_sp,0,0,C_in,d_in,param_sp_uzawa);
        end
        
        %Coordination :
        p = p + rho*sum(P);
        
        %R�cup�ration de la solution a l'it�ration k :
        P_sol = [P_sol,P];
        
        %R�cup�ration du prix � l'it�ration k:
        Prix = [Prix,p];
        
        %Incr�mentation du nombre d'it�rations:
        k = k + 1;
    end
    
    %Calcul de la valeur optimale
    J = sum(a.*(P-P0).^2 + b);
    
    t = toc(start);
end

