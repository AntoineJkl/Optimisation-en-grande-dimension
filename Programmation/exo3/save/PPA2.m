clear all

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

% contraintes
theta = [zeros(1, N); ones(1, N); -ones(1, N); -eye(N)]; %theta i en colonnes
v = [De; 1; -1; zeros(N, 1)];

%-----------------------Initialisation des algorithmes
tic

u = ones(N, 1)/N; % solution initiale

%u = [1; zeros(N-1, 1)]; % solution initiale
err = @(x, y) norm(x-y, 2)/norm(x, 2); % fonction d'erreur

% problème auxilaire
ppa = {};
ppa.k = 1;
ppa.kmax = 1000;
ppa.seuil = 1e-6;

% prédiction par les prix
prix = {};
prix.p = zeros(size(v));
prix.it = 1;
prix.it_max = 1000;
prix.seuil = 1e-5;
prix.eps = .01;

%-----------------------Resolution du problème

act = @(u) [u'*Q*u; sum(u); -sum(u); -u];

while (ppa.k<=2 || ((err(u, u_prec) + err(ppa.p, ppa.p_prec)+ norm(theta*u-v))>ppa.seuil && ppa.k<ppa.kmax))
    u_prec = u;
    ppa.p_prec = prix.p;
    
    ppa.A = eye(N);
    ppa.b = -e;
    
    theta(1,:) = u'*Q;
    
    % decomposition par les prix
    prix.it = 1;
    while(prix.it <= 1 || err(prix.p, prix.p_prec) > prix.seuil && prix.it < prix.it_max) 
        prix.p_prec = prix.p;
        
        for i = 1:N
            u(i) = (-theta(:,i)'* prix.p + ppa.b(i))/(2*ppa.A(i,i));
        end
        theta(1,:) = u'*Q;
        
        prix.p = max(0, prix.p + prix.eps * (act(u) - v));
        prix.it = prix.it + 1;
    end
    
    ppa.p = prix.p;
    ppa.k = ppa.k + 1;
end
toc

u
e'*u