function [u,v,p,k,J,u_stock] = DecompositionPrix_par(N,A,C,param)
    %N : taille de l'instance du problème
    %A : matrice de la fonction objective de taille (N+1,2,2)
    %C : matrice de contraintes associée aux sous-problèmes de taille (N,2)
    %rho : pas de gradient pour la coordination
    %eps : précision de la décomposition par les prix
    %kmax : nombre d'itérations max
    
    
    %Pour ajouter l'algorithme d'Arrow
    addpath('..\Algorithme');
    
    %Paramètres de la decomposition par prix
    if ismember('rho',fieldnames(param)) ; rho=param.rho ; else rho=1 ; end
    if ismember('eps',fieldnames(param)) ; eps=param.eps ; else eps=10^(-3) ; end
    if ismember('kmax',fieldnames(param)) ; kmax=param.kmax ; else kmax=100 ; end
    if ismember('algo',fieldnames(param)) ; algo=param.algo ; else algo='fmincon' ; end
    
    %Initialisation generale:
    k = 1;            %Iteration
    u = zeros(N+1,2); %Solution u
    v = zeros(N,2);   %Solution v
    p = zeros(1,2);   %Prix
    u_stock=zeros(kmax+2,N+1,2);
    u_stock(1,:,:)=u;
    
    %Initialisation des hyperparamètres d Arrow:
    if ismember('rho_sp1',fieldnames(param)) ; rho_sp1=param.rho_sp1 ; else rho_sp1=0.1 ; end
    if ismember('rho_sp2',fieldnames(param)) ; rho_sp2=param.rho_sp2 ; else rho_sp2=2 ; end
    if ismember('eps_sp',fieldnames(param)) ; eps_sp=param.eps_sp ; else eps_sp=10^(-10) ; end
    if ismember('kmax_sp',fieldnames(param)) ; kmax_sp=param.kmax_sp ; else kmax_sp=50000 ; end
    lambda_ini= zeros(N,1);
    mu_ini = zeros(N,5,1);
    
    tic;
    
    while( k <= 2 || ((norm(u - u_prec,2) > eps) && k <= kmax))
        %Initialisation de la valeur optimale
        J=0;
        u_prec = u;
    
        %Décomposition des N premiers sous-problèmes :
        parfor i = 1:N
            A_sp = zeros(4,4) ; %u1,u2,v1,v2
            A_sp(3,3)=(1/2)*A(i,1,1);  A_sp(4,4)=(1/2)*A(i,2,2);
            
            b_sp = zeros(4,1); %u1,u2,v1,v2
            b_sp(1:2,1)=-p;
            
            %Contrainte d inegalite
            C_in=[1 0 0 0 ; -1 0 0 0 ; 0 1 0 0 ; 0 -1 0 0 ; -1 0 -1 0];
            d_in=[C(i,1) ; 0 ; C(i,2)-C(i,1) ; 0 ; -C(i,1)];
            
            %Contrainte d inegalite
            C_eq=ones(1,4);
            d_eq=C(i,2);
            
            u_ini = [u(i,:)  v(i,:)]';
            
            if strcmp(algo,'arrow')
                %Resolution par Arrow-Hurwicz
                param_sp = struct('rho1', rho_sp1, ...
                        'rho2', rho_sp2, ...
                        'mu_ini' , mu_ini(i,:,:)' , ...
                        'lambda_ini' , lambda_ini(i,:) , ...
                        'eps', eps_sp, ...
                        'kmax', kmax_sp, ...
                        'U_ini' , u_ini);           

                [temp,lambda_ini(i,:),mu_ini(i,:,:),~] = ArrowHurwicz(A_sp,b_sp,C_eq,d_eq,C_in,d_in,param_sp);
                u(i,:)=temp(1:2);
                v(i,:)=temp(3:4); 
                
            else
                %Resolution par la méthode du point intérieur
                fun = @(var) Juv(var,A_sp,b_sp);
                options = optimoptions('fmincon','Display', 'off','GradObj','on');
                ub=[] ; lb=[] ; noncol= [] ; 
                temp = fmincon(fun,u_ini,C_in,d_in,C_eq,d_eq,lb,ub,noncol,options);
                u(i,:)=temp(1:2);
                v(i,:)=temp(3:4);
                
            end
            

            
            %Incrementation de la valeur objective J pour les N premiers
            %sous-problèmes
            J=J+temp'*A_sp*temp;
            
        end
        
        %Décomposition du N+1ème sous-problème
        u(N+1,1)=p(1,1)/A(N+1,1,1); 
        u(N+1,2)=p(1,2)/A(N+1,2,2);
        
        %Incrémentation de la valeur objective pour le N+1ème
        %sous-problèmes
        A_sp=zeros(2,2);
        A_sp(1,1)=A(N+1,1,1) ; A_sp(2,2)=A(N+1,2,2);
        J=J+(1/2)*(u(N+1,:)*A_sp*u(N+1,:)');
        
        %Coordination:
        p = p + rho.*(sum(u(1:N,:),1)-u(N+1,:));
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
        u_stock(k,:,:)=u;
    end

    toc;
end

