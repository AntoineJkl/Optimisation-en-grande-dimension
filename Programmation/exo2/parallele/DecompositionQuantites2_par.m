function [u,v,w,k,J,u_stock] = DecompositionQuantites2_par(N,A,C,param)
    %N : taille de l'instance du problème
    %A : matrice de la fonction objective de taille (N+1,2,2)
    %C : matrice de contraintes associée aux sous-problèmes de taille (N,2)
    %rho : pas de gradient pour la coordination
    %eps : précision de la décomposition par les prix
    %kmax : nombre d'itérations max
    
    tic;
       
    %Pour ajouter l'algorithme d'Uzawa 
    addpath('..\Algorithme');
    
    %Paramètres de la decomposition par quantité
    if ismember('rho',fieldnames(param)) ; rho=param.rho ; else rho=1 ; end
    if ismember('eps',fieldnames(param)) ; eps=param.eps ; else eps=10^(-3) ; end
    if ismember('kmax',fieldnames(param)) ; kmax=param.kmax ; else kmax=100 ; end
    if ismember('algo',fieldnames(param)) ; algo=param.algo ; else algo='fmincon' ; end
    
    %Initialisation generale:
    k = 1; %Iteration
    u = zeros(N+1,2); %Solution u
    v = zeros(N,2);   %Solution v
    p= zeros(N+1,2); %Multiplicateur de Lagrange
    lambda_ini= zeros(N,1);
    mu_ini = zeros(N,1);
    u_stock=zeros(kmax+2,N+1,2);
    u_stock(1,:,:)=u;
    
    %Initialisation des Allocations
    w = C;
    
    %Initialisation des hyperparamètres d'Uzawa:
    if ismember('rho_sp',fieldnames(param)) ; rho_sp=param.rho_sp ; else rho_sp=3 ; end
    if ismember('eps_sp',fieldnames(param)) ; eps_sp=param.eps_sp ; else eps_sp=10^(-10) ; end
    if ismember('kmax_sp',fieldnames(param)) ; kmax_sp=param.kmax_sp ; else kmax_sp=50000 ; end
    
    while( k <= 2 || ((norm(u - u_prec,2) > eps) && k <= kmax))
        %Initialisation de la valeur optimale
        J=0;
        
        u_prec = u;
        
        %Décomposition des N premiers sous-problèmes :
        parfor i = 1:N
            A_sp = (1/2).*reshape(A(i,:,:),2,2); %v1,v2
            
            b_sp = zeros(2,1); %v1,v2

            %Contrainte d inegalite
            C_in=[-1 0 ];
            d_in=-C(i,1)+w(i,1);
            
            %Contrainte d egalite
            C_eq=ones(1,2);
            d_eq=C(i,2)-w(i,1)-w(i,2);
            
            if strcmp(algo,'uzawa')
                %Resolution par Uzawa
                param_sp = struct('rho', rho_sp, ...
                        'mu_ini' , mu_ini(i,:) , ...
                        'lambda_ini' , lambda_ini(i,:) , ...
                        'eps', eps_sp, ...
                        'kmax', kmax_sp); 

                [v(i,:),lambda,mu,~] = Uzawa(A_sp,b_sp,C_eq,d_eq,C_in,d_in,param_sp);
                u(i,:)=w(i,:);
                p(i,:)=[ mu-lambda , -lambda];
                %p(i,1)=mu-lambda; p(i,2)=-lambda;
                lambda_ini(i,:)=lambda ; mu_ini(i,:)=mu;
           elseif strcmp(algo,'explicite')
               %Resolution par formule explicite
                if A(i,1,1)*(C(i,1)-w(i,1))+(C(i,1)+w(i,2)-C(i,2))*A(i,2,2) >= 0
                    mu=A(i,1,1)*(C(i,1)-w(i,1))+(C(i,1)+w(i,2)-C(i,2))*A(i,2,2);
                    lambda=(C(i,1)+w(i,2)-C(i,2))*A(i,2,2);
                    v(i,:) = [C(i,1)-w(i,1) , -lambda/A(i,2,2)];
                else
                    mu=0;
                    lambda=(w(i,1)+w(i,2)-C(i,2))*(A(i,1,1)*A(i,2,2))/(A(i,1,1)+A(i,2,2));
                    v(i,:)=[-lambda/A(i,1,1) , -lambda/A(i,2,2)];
                end
                u(i,:)=w(i,:);
                p(i,:)=[ mu-lambda , -lambda];
            else
                %Resolution par la méthode du point intérieur
                fun = @(var) Juv(var,A_sp,b_sp);
                options = optimoptions('fmincon','Display', 'off','GradObj','on');
                ub=[] ; lb=[] ; noncol= [] ; 
                [v(i,:),~,~,~,lagrange] = fmincon(fun,v(i,:)',C_in,d_in,C_eq,d_eq,lb,ub,noncol,options);
                lambda=lagrange.eqlin ;mu=lagrange.ineqlin ;
                p(i,:)=[ mu-lambda , -lambda];
                u(i,:)=w(i,:);
            end
            
            %Incrementation de la valeur objective J pour les N premiers
            %sous-problèmes
            J=J+v(i,:)*A_sp*v(i,:)';
            
        end
        
        %Décomposition du N+1ème sous-problème
        u(N+1,:)=sum(w,1);
        p(N+1,1)=u(N+1,1)*A(N+1,1,1) ; p(N+1,2)=u(N+1,2)*A(N+1,2,2);
        
        %Incrémentation de la valeur objective pour le N+1ème
        %sous-problèmes
        A_sp=zeros(2,2);
        A_sp(1,1)=A(N+1,1,1) ; A_sp(2,2)=A(N+1,2,2);
        J=J+(1/2)*(u(N+1,:)*A_sp*u(N+1,:)');
        
        %Coordination:
        w = w + rho*(p(1:N,:)-repmat(p(N+1,:),N,1));
        w(:,1)=min(max(0,w(:,1)),C(:,1));
        w(:,2)=min(max(0,w(:,2)),C(:,2)-C(:,1));
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
        u_stock(k,:,:)=u;
        
    end
    toc;


end

