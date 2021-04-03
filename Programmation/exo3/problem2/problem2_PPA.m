function [ u ] = problem2_PPA(Q, e, De, param)
%PROBLEM2_PPA Resolution du probleme 2 par la methode du probleme auxiliaire
% Q: matrice de covariance
% e: esperance
% De: dose de risque maximale
% param: struct avec | alpha : force convexite auxiliaire
%                    | eps : coefficient auxiliaire 
%                    | beta : actualisation des prix 
%                    | kmax : nombre d'it maximal
%                    | seuil :  seuil d'erreur
    
    %----------------------- Verification que Q est SDP
    N = length(e);
    if (min(eig(Q))<=0)
        disp(['Attention ! Matrice Q non SDP : l_min = ' num2str(min(eig(Q)))]);
        return;
    end
    
    %----------------------- Parametres de l'algorithme d'Uzawa
    addpath('../../Algorithme')
    param_Uzawa = struct('rho', .1, ...
                         'mu', 0, ...
                         'lambda', 0);
    mu = zeros(N,1); % vecteur des multiplicateurs
    
    %----------------------- Variables du probleme
    QA = diag(diag(Q));
    QC = Q - QA;

    u = ones(N, 1)/N; %solution initiale
    U = zeros(N, param.kmax); %solutions a chaque instant
    
    err = @(x, y) norm(x-y, 2)/norm(x, 2); % fonction d'erreur

    % parametres de la decomposition par les prix
    prix = struct('p', [0 0], ...      %prix
                  'p_prec', [0 0], ... %prix a l'iteration precedente
                  'beta', param.beta); %terme de convexite

    % parametres du probleme auxilaire
    ppa = struct('k',1, ...
                 'kmax', param.kmax, ...
                 'seuil', 1e-5, ...
                 'eps', param.eps, ...
                 'alpha', param.alpha);

    %-----------------------Resolution du probleme
    k = 1;
    while (k<=2 || ((err(u, u_prec) + err(prix.p, prix.p_prec))>ppa.seuil && k<ppa.kmax))
        % sauvegarde de l'iteration precedente
        U(:, k) = u; 
        prix.p_prec = prix.p;
        u_prec = u;

        % actualisation des termes du probleme
        ppa.A = ppa.alpha(k) * eye(N)/2 + ppa.eps(k) * prix.p(1) * QA;
        ppa.b = ppa.eps(k)* e + ppa.alpha(k) * u - ppa.eps(k) * prix.p(1)*2*QC*u - ppa.eps(k) * prix.p(2) * ones(N, 1);

        % resolution des sous-problemes par Uzawa
        for i = 1:N
            param_Uzawa.Mu = mu(i); 
            [u(i), ~, mu(i)] = Uzawa(ppa.A(i,i), ppa.b(i), 0, 0, -1, 0, param_Uzawa);
        end
        
        % coordination des prix
        prix.p(1) = max(0, prix.p(1) + prix.beta * ppa.eps(k) * (u'*Q*u - De));
        prix.p(2) = prix.p(2) + prix.beta * ppa.eps(k) * (sum(u) - 1);

        k = k + 1;
    end
    
    u = U(:, 1:k-1);
end

