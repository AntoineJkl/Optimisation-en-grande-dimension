clear 
%----------------------- Données
% jouet 1
N = 2; % nombre d'actions
Q = [2 1; 1 2];
e = [4; 5];
De = 3/2;
sol_exacte = [1; 1]/2;
% Minimize[{4x + 5y, x+y=1,x>=0,y>=0, 2*x^2 + 2*y^2 + 1*x*y <= 3}, {x, y}]

% jouet 2
% N = 2; % nombre d'actions
% Q = [3/2 0; 0 2];
% e = [1; 2];
% De = 1;
% sol_exacte = ([4; 3] + sqrt(2)*[-1; 1])/7;
% Minimize[{- 1x - 2y, x+y = 1, x>=0, y>=0, 3/2*x^2 + 2*y^2 + 0*x*y <= 1}, {x, y}]

% ------------ Résolution du problème
param = struct('alpha', @(k) .1, ... #convexite auxiliaire
               'eps', @(k) .1, ... #auxiliaire 
               'beta', 10, ...  #actualisation des prix 
               'kmax', 1000,...
               'seuil', 1e-6);
tic
u = problem_2(Q, e, De, param);
toc

verif_sol = @(x) struct('J', [e'*x, e'*sol_exacte], ...
                        'l', [sum(x) 1], ...
                        'De',[x'*Q*x, De], ...
                        'u', [x', sol_exacte']);

verif_sol(u(:,end))

%% Affichage de la solution
figure(1);
color = ['m' 'c'];
for i = 1:length(color)
    plot(u(i, :), 'Color', color(i));
    hold on;
    plot(xlim, sol_exacte(i)*[1 1], '--','Color', color(i), 'LineWidth', .2);
end
hold off
legend('u_1', 'u_1 exact', 'u_2', 'u_2 exact');
title('Visualisation de la convergence');
ylim([0 1]);