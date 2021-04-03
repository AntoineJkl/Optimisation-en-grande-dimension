function [ u ] = problem_1(Q, e, Re, param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    addpath('../Algorithme')

    N = length(e);
    if (min(eig(Q))<=0)
        disp(['Attention ! Matrice Q non SDP : l_min = ' num2str(min(eig(Q)))]);
        return;
    end
    
    param_Uzawa = struct('rho', .1, ...
                         'Lambda', 0, ...
                         'Mu', 0);
    mu = zeros(N,1);
    
    kmax = param.kmax;
    %----------------------- Variables du problème
    % décomposition de l'objectif
    QA = diag(diag(Q));
    QC = Q - QA;

    U = zeros(N, kmax);
    u = ones(N, 1)/N; % solution initiale
    err = @(x, y) norm(x-y, 2); % fonction d'erreur

    % décomposition par les prix
    prix = struct('p', [0 0], ...      %prix
                  'p_prec', [0 0], ... %prix a l'iteration precedente
                  'beta', param.beta);         %terme de convexite

    % problème auxilaire
    ppa = struct('k',1, ...
                 'kmax', param.kmax, ...
                 'seuil', param.seuil, ...
                 'eps', param.eps, ...   %terme actualisation
                 'alpha', param.alpha);

    %-----------------------Resolution du problème
    k = 1;
    while (k<=2 || ((err(u, u_prec) + err(prix.p, prix.p_prec))>ppa.seuil && k<ppa.kmax))
        U(:, k) = u;
        prix.p_prec = prix.p;
        u_prec = u;

        ppa.A = ppa.alpha(k) * eye(N)/2 + ppa.eps(k) * QA/2;
        ppa.b = ppa.alpha(k) * u - ppa.eps(k)*QC*u + ppa.eps(k)*prix.p(1)*e - ppa.eps(k) * prix.p(2) * ones(N, 1);

        % decomposition par les prix
        
        for i = 1:N
            param_Uzawa.Mu = mu(i);
            [u(i), ~, mu(i)] = Uzawa(ppa.A(i,i), ppa.b(i), 0, 0, -1, 0, param_Uzawa);
        end

        prix.p(1) = max(0, prix.p(1) + prix.beta * ppa.eps(k) * (-e'*u + Re));
        prix.p(2) = prix.p(2) + prix.beta * ppa.eps(k) * (sum(u) - 1);
        
        k = k + 1;
    end
    u = U(:, 1:k-1);
end

