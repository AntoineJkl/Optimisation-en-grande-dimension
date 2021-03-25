function [u,v,w,k,J] = DecompositionQuantites2(N,A,C,rho,eps,kmax)
    %N : taille de l'instance du problème
    %A : matrice de la fonction objective de taille (N+1,2,2)
    %C : matrice de contraintes associée aux sous-problèmes de taille (N,2)
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
    w = ones(N,2)./(2*N); %Quantités
    p= zeros(N+1,2);
    
    %Initialisation des hyperparamètres de Uzawa ou Arrow:
    rho_sp = 0.1;
    eps_sp = 10^(-10);
    kmax_sp = 50000;
    
    while( k <= 2 || ((norm(u - u_prec,2)/norm(u,2) + norm(w - w_prec,2)/norm(w,2) + norm(v - v_prec,2)/norm(v,2) > eps) && k <= kmax))
        %Initialisation de la valeur optimale
        J=0;
        
        u_prec = u;
        v_prec = v;
        w_prec = w;
        
        %Décomposition des N premiers sous-problèmes :
        for i = 1:N
            A_sp = (1/2).*reshape(A(i,:,:),2,2); %v1,v2
            
            b_sp = zeros(2,1); %v1,v2

            C_in=[-1 0 ];
            d_in=-C(i,1)+w(i,1);
            
            C_eq=ones(1,2);
            d_eq=C(i,2)-w(i,1)-w(i,2);
            
            mu_ini=0;
            lambda_ini=0;
            
            [v(i,:),lambda,mu,~] = Uzawa(A_sp,b_sp,C_eq,d_eq,C_in,d_in,rho_sp,mu_ini,lambda_ini,eps_sp,kmax_sp);
            u(i,:)=w(i,:);
            p(i,1)=mu-lambda; p(i,2)=-lambda;
            
            %Incrementation de la valeur objective J pour les N premiers
            %sous-problèmes
            J=J+v(i,:)*A_sp*v(i,:)';
            
        end
        
        %Décomposition du N+1ème sous-problème
        u(N+1,:)=-sum(w,1);
        p(N+1,1)=u(N+1,1)*A(N+1,1,1) ; p(N+1,2)=u(N+1,2)*A(N+1,2,2);
        
        %Incrémentation de la valeur objective pour le N+1ème
        %sous-problèmes
        A_sp=zeros(2,2);
        A_sp(1,1)=A(N+1,1,1) ; A_sp(2,2)=A(N+1,2,2);
        J=J+(1/2)*(u(N+1,:)*A_sp*u(N+1,:)');
        
        %Coordination:
        for i = 1:(N)
            w(i,:) = w(i,:) + rho.*(p(i,:)-(1/(N+1)).*sum(p,1));
            w(i,1)=min(max(0,w(i,1)),C(i,1));
            w(i,2)=min(max(0,w(i,2)),C(i,2)-C(i,1));
        end
        %Incrementation du nombre d'iterations:
        k = k + 1;
    end

    for i=1:N
        if u(i,1) > C(i,1) || u(i,2) > C(i,2)-C(i,1) || u(i,1)< 0 || u(i,2) <0
            disp('Contrainte sur l ensemble admissible de u non respecté'); 
        end
    end
    toc;


end

