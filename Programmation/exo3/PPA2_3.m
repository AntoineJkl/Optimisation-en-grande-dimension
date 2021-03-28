clear all
addpath('../Algorithme')
%----------------------- Données
% jouet 1
N = 2; % nombre d'actions
Q = [2 1; 1 2];
e = [4; 5];
De = 3/2;

% jouet 2
N = 2; % nombre d'actions
Q = [5/3 3; 3 10/3];
e = [5/4; 3];
De = 3;
% [1 2]/3 -> 29/12

%----------------------- Variables du problème

%-----------------------Initialisation des algorithmes
tic

QA = diag(diag(Q));
QC = Q - QA;

u = ones(N, 1)/N; % solution initiale

err = @(x, y) norm(x-y, 2)/norm(x, 2); % fonction d'erreur

% problème auxilaire
ppa = struct('k',1, ...
             'kmax', 1000, ...
             'seuil', 1e-6);

% prédiction par les prix
prix = struct('p', [0 0 0], ... %zeros(size(v) + [1 0]), ...
              'it', 1, ...
              'it_max', 1000, ...
              'seuil', 1e-6, ...
              'eps', .1);

%-----------------------Resolution du problème

act = @(u) [u'*Q*u; sum(u); -sum(u); -u];

eps = @(k) .1;
k = 1;
while (k<=2 || ((err(u, u_prec) + err(prix.p, prix.p_prec))>ppa.seuil && k<ppa.kmax))
    prix.p_prec = prix.p;
    u_prec = u;
    prix.it = 1;
    
    ppa.A = eye(N)/2 + eps(k) * prix.p(1) * QA;
    ppa.b = eps(k)* e + u - eps(k) * prix.p(1)*2*QC*u - eps(k) * prix.p(2) * ones(N, 1) + eps(k) * prix.p(3) * ones(N, 1);
    
    % decomposition par les prix
    %u = Uzawa(ppa.A, ppa.b, ones(1,N), 1, -eye(N), zeros(N, 1), .1);
    for i = 1:N
        u(i) = Uzawa(ppa.A(i,i), ppa.b(i), 0, 0, -1, 0, 1);
    end
    
    %u(N) = sum(u(2:N-1));
    prix.p(1) = max(0, prix.p(1) + eps(k) * prix.eps * (u'*Q*u - De));
    prix.p(2) = max(0, prix.p(2) + eps(k) * prix.eps * (sum(u) - 1));
    prix.p(3) = max(0, prix.p(3) + eps(k) * prix.eps * (-sum(u) + 1));
    
    %prix.p(2:end) = max(0, prix.p(2:end) + eps(k) * prix.eps * (theta * u - v));
    k = k + 1;
end
toc

u
e'*u