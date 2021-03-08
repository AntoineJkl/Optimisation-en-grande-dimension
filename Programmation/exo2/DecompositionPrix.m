function [u,p,k,J] = DecompositionPrix(N,A,C,rho,eps,kmax)
    %N : taille de l'instance du problème
    %A : matrice de la fonction objective de taille (N+1,2,2)
    %C : matrice de contraintes associée aux sous-problèmes de taille (N+1,2)
    %rho : pas de gradient pour la coordination
    %eps : précision de la décomposition par les prix
    %kmax : nombre d'itérations max
    
    tic;
       
    %Pour ajouter les algorithmes (Uzawa et Arrow)
    addpath('..\Algorithme');
    
    %Initialisation generale:
    k = 1; %Iteration
    u = zeros(N+1,2); %Solution u
    v = zeros(N,2);   %Solution v
    p = zeros(1,2); %Prix
    
    %Initialisation des hyperparamètres de Uzawa ou Arrow:
    rho_sp = 0.1;
    eps_sp = 10^(-5);
    kmax_sp = 10000;
    
    while( k <= 2 || ((norm(u - u_prec,2)/norm(u,2) + norm(p - p_prec,2)/norm(p,2) + norm(v - v_prec,2)/norm(v,2) > eps) && k <= kmax))
        u_prec = u;
        v_prec = v;
        p_prec = p;
        
        %Décomposition des N premiers sous-problèmes :
        for i = 1:N
            A_sp = zeros(4,4) ; %u1,u2,v1,v2
            A_sp(3,3)=(1/2).*A(i,1,1);  A_sp(4,4)=(1/2).*A(i,2,2);
            
            b_sp = zeros(4,1); %u1,u2,v1,v2
            b_sp(1:2,1)=p;
            
            C_in=[1 0 0 0 ; -1 0 0 0 ; 0 1 0 0 ; 0 -1 0 0 ];
            d_in=[C(i,1) ; 0 ; C(i,2)-C(i,1) ; 0 ];
            
            C_eq=ones(1,4);
            d_eq=C(i,2);
            
            mu_ini=zeros(4,1);
            lambda_ini=0;
            
            [temp,~,~,~] = Uzawa(A_sp,b_sp,C_eq,d_eq,C_in,d_in,rho_sp,mu_ini,lambda_ini,eps_sp,kmax_sp);
            u(i,:)=temp(1:2);
            v(i,:)=temp(3:4);
            
        end
        
        %Décomposition du N+1ème sous-problème
        u(N+1,:)=p/C(N+1,:);
        
        %Coordination:
        p = p + rho.*(sum(u(1:N,:),1)-u(N+1,:));
        
        %Incrementation du nombre d'iterations:
        k = k + 1;
    end
    
    %Calcul de la valeur optimale
    J = 0; %1/2*u'*A*u - b'*u

    toc;
end

