function [u,p,k,J,t2,U_final,Mu_final] = DecompositionPrix(N,A,b,C,parametres,parametres_sousproblemes)
    t1=tic;
    
    %Pour ajouter les algorithmes (Uzawa et Arrow)
    addpath('..\Algorithme');
    
    %Paramètres;
    if ismember('rho',fieldnames(parametres)) ; rho=parametres.rho ; else rho=0.05 ; end
    if ismember('eps',fieldnames(parametres)) ; eps=parametres.eps ; else eps=10e-5 ; end
    if ismember('kmax',fieldnames(parametres)) ; kmax=parametres.kmax ; else kmax=5000 ; end
    if ismember('PrintIt',fieldnames(parametres)) ; PrintIt=parametres.PrintIt ; else PrintIt=false ; end
    if ismember('bigU',fieldnames(parametres)) ; bigU=parametres.bigU ; else bigU=false ; end
    if ismember('bigMu',fieldnames(parametres)) ; bigMu=parametres.bigU ; else bigMu=false ; end
    if ismember('rho_sp',fieldnames(parametres_sousproblemes)) ; rho_sp=parametres_sousproblemes.rho_sp_uzawa ; else rho_sp=0.01 ; end
    if ismember('eps_sp',fieldnames(parametres_sousproblemes)) ; eps_sp=parametres_sousproblemes.eps_sp ; else eps_sp = 10e-4 ; end
    if ismember('kmax_sp',fieldnames(parametres_sousproblemes)) ; kmax_sp=parametres_sousproblemes.kmax_sp ; else kmax_sp = 3000 ; end
    
    
    %Initialisation generale:
    U_final=[]; %Vecteur des solutions à chaque iteration i
    Mu_final=[]; %Vecteur des multiplicateurs à chaque iteration i
    k = 1; %Iteration
    u = zeros(N,1); %Solution
    p = zeros(N,1); %Prix
    critere = 0; %Critere d'arret
    
    while( k <= 2 || (~critere && k <= kmax))
        
        %Affichage de l'iteration courante:
        if PrintIt
            disp(['Iteration - Prix: ',num2str(k)]);
        end
        
        %Solution et prix de l'etape precedente:
        u_prec = u;
        p_prec = p;
        
        %Décomposition:
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
        
        %Mise à jour du critere d'arret:
        critere = (norm(u - u_prec,2) < eps);
        
        if bigU
        %Mise à jour des solutions calculées
            U_final=[U_final,u];
        end
        
        if bigMu
        %Mise à jour des multiplicateurs calculées
            Mu_final=[Mu_final,p];
        end
    end
    
    %Calcul de la valeur optimale
    J = 1/2*u'*A*u - b'*u;
    
    t2=toc(t1);
end

