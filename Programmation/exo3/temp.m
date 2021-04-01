% Fichier de test pour des instances de taille N 
%% Problème 1

clear all
count = 1;
for N = 3:3;
    Q = 2 * eye(N) + diag(ones(N-1, 1), -1) + diag(ones(N-1, 1), 1);
    e = 2 + (1:N)'/2;
    Re = sum(e);
    
    count = count * (min(eig(Q)) > 0);
    
    param = struct('alpha', @(k) 4, ... #convexite auxiliaire
                   'eps', @(k) .3, ... #auxiliaire 
                   'beta', 1, ...  #actualisation des prix 
                   'kmax', 1000,...
                   'seuil', 1e-6);

    u = problem_1(Q, e, Re, param);
    
end
rep = u(:, end);
[e'*rep Re]
sum(rep)

%% Problème 2
for N = 2:2;
    Q = 2 * eye(N) + diag(ones(N-1, 1), -1) + diag(ones(N-1, 1), 1);
    e = 2 + (1:N)'/2;
    De = sum(sum(Q))/2;
        
    param = struct('alpha', @(k) 1, ... #convexite auxiliaire
                   'eps', @(k) .1, ... #auxiliaire 
                   'beta', 500, ...  #actualisation des prix 
                   'kmax', 1000,...
                   'seuil', 1e-6);

    u = problem_2(Q, e, De, param);
    
end
rep = u(:, end);
[e'*rep De]
sum(rep)
