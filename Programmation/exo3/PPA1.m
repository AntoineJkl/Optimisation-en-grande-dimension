clear all

%----------------------- Données
% jouet 1
N = 2; % nombre d'actions
Q = [2 1; 1 2];
e = [4; 5];
Re = 9/2;

% jouet 2
N = 2; % nombre d'actions
Q = [3/5 3; 3 10/3];
e = [5/4; 3];
Re = 3/2;

%----------------------- Variables du problème
% décomposition de l'objectif
QA = diag(diag(Q)); %partie additive
QC = Q - QA; %partie couplante

% contraintes
theta = [-e'; ones(size(e')); -ones(size(e'))]; %theta i en colonnes
v = [-Re; 1; -1];

%-----------------------Initialisation des algorithmes
tic

u = ones(N, 1)/N; % solution initiale
err = @(x, y) norm(x-y, 2)/norm(x, 2); % fonction d'erreur

% problème auxilaire
ppa = {};
ppa.k = 1;
ppa.kmax = 1000;
ppa.seuil = 1e-6;

% prédiction par les prix
prix = {};
prix.p = [0; 0; 0];
prix.it = 1;
prix.it_max = 500;
prix.seuil = 1e-6;
prix.eps = .1;

%-----------------------Resolution du problème
while (ppa.k<=2 || ((err(u, u_prec) + err(ppa.p, ppa.p_prec))>ppa.seuil && ppa.k<ppa.kmax))
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
    ppa.k = ppa.k + 1;
end
toc

u