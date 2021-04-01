function [ u ] = problem_1_dorine(Q, e, Re)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    N = length(e);
    if (min(eig(Q))<=0)
        disp(['Attention ! Matrice Q non SDP : l_min = ' num2str(min(eig(Q)))]);
        return;
    end
    
    param_Uzawa = struct('rho', .1);
    
    %----------------------- Variables du problème
    % décomposition de l'objectif
    QA = diag(diag(Q)); %partie additive
    QC = Q - QA; %partie couplante

    % contraintes
    theta = [-e'; ones(size(e')); -ones(size(e')); -eye(N)]; %theta i en colonnes
    v = [-Re; 1; -1; zeros(N,1)];
    
    kmax = 1000;
    U = zeros(N, kmax); % solution initiale 
   
    u = ones(N, 1)/N; % solution initiale
    
    %-----------------------Initialisation des algorithmes
    % problème auxilaire
    ppa = struct('kmax', kmax, ...
                 'seuil', 1e-6);
    
    % prédiction par les prix
    prix = struct('p', zeros(size(v)), ...
                  'it', 1, ...
                  'it_max', 500, ...
                  'seuil', 1e-6, ...
                  'eps', .1);

    err = @(x, y) norm(x-y, 2); % fonction d'erreur
    
    k = 1;
    %-----------------------Resolution du problème
    while (k<=2 || ((err(u, u_prec) + err(ppa.p, ppa.p_prec))>ppa.seuil && k<ppa.kmax))
        U(:, k) = u;
        u_prec = u;
        ppa.p_prec = prix.p;

        ppa.A = QA;
        ppa.b = -1/2 * Q * u + QA * u;

        % prédiction par les prix
        prix.it = 1;
        while(prix.it <= 1 || err(prix.p, prix.p_prec) > prix.seuil && prix.it < prix.it_max) 
            prix.p_prec = prix.p;

            for i = 1:N
                u(i) = (-theta(:,i)'* prix.p + ppa.b(i))/(2*ppa.A(i,i));
            end

            prix.p = max(0, prix.p + prix.eps * (theta * u - v));
            prix.it = prix.it + 1;
        end

        ppa.p = prix.p;
        k = k + 1;
    end
    u = U(:, 1:k-1);
end

