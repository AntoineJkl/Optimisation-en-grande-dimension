function [U,Lambda,Mu,k] = Uzawa(A,b,C_eq,d_eq,C_in,d_in,param)
%UZAWA Summary of this function goes here
%   Detailed explanation goes here

    %Arguments d'entree
    if ismember('rho',fieldnames(param)) ; rho=param.rho ; else rho=0.1 ; end
    if ismember('mu_ini',fieldnames(param)) ; mu_ini=param.mu_ini ; else mu_ini=zeros(size(C_in,1),1) ; end
    if ismember('lambda_ini',fieldnames(param)) ; lambda_ini=param.lambda_ini ; else lambda_ini=zeros(size(C_eq,1),1) ; end
    if ismember('eps',fieldnames(param)) ; eps=param.eps ; else eps=10^(-5) ; end
    if ismember('kmax',fieldnames(param)) ; kmax=param.kmax ; else kmax=10000 ; end
    
    %Initialisation
    k=1;
    U=[];
    Mu=mu_ini;
    Lambda=lambda_ini;
    
    %Résolution problème
    while ( k<=2 || ( norm(U-prec,2)  > eps ) && ( k<=kmax ) )  

        %Stockage de u précédent
        prec=U;
        
        %Annulation du gradient du lagrangien
        U = (2*A)\(b-C_in'*Mu-C_eq'*Lambda); 
        
        %Mis à jour des multiplicateurs de Lagrange
        Mu=max(0, Mu+rho.*(C_in*U-d_in));
        Lambda=Lambda+rho.*(C_eq*U-d_eq);
        
        %Incrémentation du nombre d'itérations
        k=k+1; 
    end
    
end

