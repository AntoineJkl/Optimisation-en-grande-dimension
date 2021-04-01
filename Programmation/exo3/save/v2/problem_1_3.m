function [ u ] = problem_1_3(Q, e, Re, param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    addpath('../Algorithme')

    N = length(e);
    if (min(eig(Q))<=0)
        disp(['Attention ! Matrice Q non SDP : l_min = ' num2str(min(eig(Q)))]);
        return;
    end
    
    param_Uzawa = struct('rho', .01, ...
                         'Lambda', 0, ...
                         'Mu', 0);
    mu = zeros(N,1);
    
    
    kmax = 1000;
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
    lambda = zeros(N, 1);
    while (k<=2 || ((err(u, u_prec) + err(prix.p, prix.p_prec))>ppa.seuil && k<ppa.kmax))
        U(:, k) = u;
        prix.p_prec = prix.p;
        u_prec = u;

        ppa.A = QA; %ppa.alpha(k) * eye(N)/2 + ppa.eps(k) * QA/2;
        ppa.b = -1/2 * Q * u + QA * u; %ppa.alpha(k) * u - ppa.eps(k)*QC*u ;% + ppa.eps(k)*prix.p(1)*e - ppa.eps(k) * prix.p(2) * ones(N, 1);

        % decomposition par les prix
        it = 1;
        while (it <= 2 || (err(prix.p, p_prec) + err(lambda, lambda_prec)> 1e-6 && it < 10000))
            p_prec = prix.p;
            lambda_prec = lambda;
            theta = [e'; ones(size(e')); -eye(N)]; %theta i en colonnes
            v = [Re; 1; zeros(N,1)];
            
            for i = 1:N
                %u(i) = (lambda(i) - ppa.eps(k) * prix.p(1)*e(i) - ppa.eps(k)* prix.p(2) + - ppa.eps(k)* prix.p(3) + ppa.b(i))/(2*ppa.A(i,i));
                u(i) = (-theta(:,i)'* [prix.p'; lambda] + ppa.b(i))/(2*ppa.A(i,i));
            end

            prix.p(1) = max(0, prix.p(1) + prix.beta * ppa.eps(k) * (e'*u - Re));
            prix.p(2) = prix.p(2) + prix.beta * ppa.eps(k) * (sum(u) - 1);
            lambda = max(0, lambda + prix.beta * ppa.eps(k) * (-u));
            it = it + 1;
        end
        k = k + 1;
    end
    u = U(:, 1:k-1);
end

