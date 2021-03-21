addpath('..\Algorithme');

%----------------------- Données
% jouet 1
N = 2; % nombre d'actions
Q = [2 1; 1 2];
e = [4; 5];
Re = 9/2;

% % jouet 2
% N = 2; % nombre d'actions
% Q = [3/5 3; 3 10/3];
% e = [5/4; 3];
% Re = 3/2;


%----------------------- Variables du problème
%décomposition de l'objectif
QA = diag(diag(Q));
QC = Q - QA;

%contraintes
theta = [-e'; ones(size(e'))]; %theta i en colonnes
v = [-Re; 1]; 

%terme auxiliaire
alpha = @(k) 1/sqrt(k); %1/2+1.2e-1;

%-----------------------
p = [0;0];
rho = 1;

eps = @(k) 1/sqrt(k);


u_prec = zeros(N, 1);

%u = [1; zeros(N-1, 1)];
%u = [6/7; 1/7];
u = [1/2; 1/2];

k = 1;
kmax = 200;
err = @(x, y) norm(x-y, 2)/norm(x, 2);
seuil = .01;

while (k<=2 || ((err(u, u_prec) + err(p, p_prec))>seuil && k<kmax))
    u_prec = u;
    
    A = alpha(k)* eye(N)/2 + eps(k) * QA/2;
    b = eps(k) * QC * u + eps(k) * theta' * p - alpha(k)*u;
    
    parfor i = 1:N
        u(i) = Uzawa(A(i,i) , b(i), 0, 0, -1, 0, rho);
    end
    
    p_prec = p;
    p(1) = max(0, p(1) + eps(k)*alpha(k) * (-e'*u + Re));
    p(2) = p(2) + eps(k)*alpha(k) * (sum(u) - 1);
    
    %beta = .4;
    %p = (1-beta) *p_prec + beta *p;
    
    k = k + 1;
end

u