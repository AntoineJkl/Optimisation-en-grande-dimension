%% Résolution du problème 2 pour une grand instance
clear all;

%-----------------------Création des données
N = 20;
Q = diag(4 + (1:N)/2) + diag((N-1):-1:1, 1)/16 + + diag(ones(N-2, 1), 2)/4;
Q = Q + triu(Q, 1)';

e = (1:N)';
De = ones(1,N)*Q*ones(N,1)/N^2 +1;

% Verification que Q est SDP
if (min(eig(Q))<=0)
    disp(['Attention ! Matrice Q non SDP : l_min = ' num2str(min(eig(Q)))]);
    break; 
end

%----------------------- Résolution du problème
sol_exacte = problem2_Solver(e, Q, De);

param = struct('alpha', @(k) 2, ... #convexite auxiliaire
               'eps', @(k) .1, ...  #auxiliaire 
               'beta', 10, ...      #actualisation des prix 
               'kmax', 1000,...
               'seuil', 1e-6);

tic
u = problem2_PPA(Q, e, De, param);
toc

%----------------------- Comparaison des solutions
str = @(x) num2str(x);
verif_sol = @(x) struct('J', e'*x , ...
                        'J_exa', e'*sol_exacte,...
                        'l', [str(sum(x)) ' == ' str(1)], ...
                        'De',[str(x'*Q*x) ' <= ' str(De)], ...
                        'u', x', ...
                        'u_exa', sol_exacte');

verif_sol(u(:,end))

%% Affichage de la solution
figure(2);
col = @(i) [i/N 0 1-i/N];

for i = 1:N
    plot(u(i, :), 'Color', col(i));
    hold on;
    plot([0 size(u,2)], sol_exacte(i)*[1 1], ':', 'Color', col(i)/2, 'LineWidth',.8);
end
hold off;

title('Evolution des u_i en fonction de l`itération');
legend('Solution PPA', 'Solution exacte');
xlabel('itération'); ylabel('u_i');
ylim([-.01 .2]);
