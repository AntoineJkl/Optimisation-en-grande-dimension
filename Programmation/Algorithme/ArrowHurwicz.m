function [U,Lambda,Mu,k] = ArrowHurwicz(A,b,C_eq,d_eq,C_in,d_in,param)
    
    %Arguments d'entree
    if ismember('rho1',fieldnames(param)) ; rho1=param.rho1 ; else rho1=0.1 ; end
    if ismember('rho2',fieldnames(param)) ; rho2=param.rho2 ; else rho2=0.1 ; end
    if ismember('mu_ini',fieldnames(param)) ; mu_ini=param.mu_ini ; else mu_ini=zeros(size(C_in,1),1) ; end
    if ismember('lambda_ini',fieldnames(param)) ; lambda_ini=param.lambda_ini ; else lambda_ini=zeros(size(C_eq,1),1) ; end
    if ismember('eps',fieldnames(param)) ; eps=param.eps ; else eps=10^(-5) ; end
    if ismember('kmax',fieldnames(param)) ; kmax=param.kmax ; else kmax=10000 ; end
    if ismember('U_ini',fieldnames(param)) ; U_ini=param.U_ini ; else U_ini=zeros(size(A,1),1); end
    
    %Initialisation
    k = 1;
    U = U_ini;
    Mu = mu_ini;
    Lambda = lambda_ini;
    
    %Resolution probleme
    while(k <= 2 || (norm( U - prec,2) > eps) && (k <= kmax) )
        %Stockage de U precedent
        prec = U;
        %Récupération de U (Formule de projection)
        U = U - rho1*(2*A*U - b + C_eq'*Lambda + C_in'*Mu);
        %Mise a jour des multiplicateurs de Lagrange
        Lambda = Lambda + rho2*(C_eq*U - d_eq);
        Mu = max(0,Mu + rho2*(C_in*U - d_in));
        %Incrementation du nombre d'iterations
        k = k + 1;
    end
    
end

