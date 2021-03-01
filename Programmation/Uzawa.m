function [U,Lambda,Mu,k] = Uzawa(A,b,d_eq,C_eq,C_in,d_in,pho,mu_ini,lambda_ini,eps,kmax)
%UZAWA Summary of this function goes here
%   Detailed explanation goes here

    %Arguments d'entrée
    switch nargin
        case 7
            mu_ini=zeros(size(C_in,1),1);
            lambda_ini=zeros(size(C_eq,1),1);
            kmax=10000;
            eps = 10^(-5);
        case 9
            kmax=10000;
            eps = 10^(-5);
        case 10
            kmax=10000;
    end
    
    %Initialisation
    k=1;
    U=[];
    Mu=mu_ini;
    Lambda=lambda_ini;
    
    %Résolution problème
    while ( k<=2 || ( norm(U-prec,2)/norm(U,2)  > eps ) && ( k<=kmax ) )  
        
        %Stockage de u précédent
        prec=U;
        
        %Annulation du gradient du lagrangien
        U= (2*A)\(b-C_in'*Mu-C_eq'*Lambda) ; 
        
        %Mis à jour des multiplicateurs de Lagrange
        Mu=max(0, Mu+pho.*(C_in*U-d_in));
        Lambda=Lambda+pho.*(C_eq*U-d_eq);
        
        %Incrémentation du nombre d'itérations
        k=k+1; 
    end
    
end

