clear all
addpath('../Algorithme')
%----------------------- Données
% jouet 1
N = 2; % nombre d'actions
Q = [2 1; 1 2];
e = [4; 5];
De = 3/2;
% [1 1]/2
% Minimize[{4x + 5y, x+y=1,x>=0,y>=0, 2*x^2 + 2*y^2 + 1*x*y <= 3}, {x, y}]

% jouet 2
% N = 2; % nombre d'actions
% Q = [3/2 0; 0 2];
% e = [1; 2];
% De = 1;
% ([4 3] + sqrt(2)*[-1 1])/7
% Minimize[{- 1x - 2y, x+y = 1, x>=0, y>=0, 3/2*x^2 + 2*y^2 + 0*x*y <= 1}, {x, y}]

%----------------------- Variables du problème

%-----------------------Initialisation des algorithmes

if (min(eig(Q))<=0)
    disp(['Attention ! Matrice Q non SDP : l_min = ' num2str(min(eig(Q)))]);
    break;
end

verif = @(u) [u'*Q*u, De;...
              sum(u), 1; ...
              u, zeros(size(u))];

tic

QA = diag(diag(Q));
QC = Q - QA;

u = ones(N, 1)/N; % solution initiale
err = @(x, y) norm(x-y, 2)/norm(x, 2); % fonction d'erreur

% décomposition par les prix
prix = struct('p', [0 0], ...      %prix
              'p_prec', [0 0], ... %prix a l'iteration precedente
              'eps', @(k) .1, ...   %terme actualisation
              'beta', 10);         %terme de convexite

% problème auxilaire
ppa = struct('k',1, ...
             'kmax', 1000, ...
             'seuil', 1e-5, ...
             'alpha', @(k) .1);

%-----------------------Resolution du problème
k = 1;

while (k<=2 || ((err(u, u_prec) + err(prix.p, prix.p_prec))>ppa.seuil && k<ppa.kmax))
    prix.p_prec = prix.p;
    u_prec = u;
    
    ppa.A = ppa.alpha(k) * eye(N)/2 + prix.eps(k) * prix.p(1) * QA;
    ppa.b = prix.eps(k)* e + ppa.alpha(k) * u - prix.eps(k) * prix.p(1)*2*QC*u - prix.eps(k) * prix.p(2) * ones(N, 1);
    
    % decomposition par les prix
    for i = 1:N
        u(i) = Uzawa(ppa.A(i,i), ppa.b(i), 0, 0, -1, 0, 1);
    end
    
    prix.p(1) = max(0, prix.p(1) + prix.beta * prix.eps(k) * (u'*Q*u - De));
    prix.p(2) = prix.p(2) + prix.beta * prix.eps(k) * (sum(u) - 1);
    
    k = k + 1;
end
toc

% Affichage de la solution
disp(['Nombre iterations : ' num2str(k)]);

verif(u)
e'*u
