clear 
% Synthaxe Wolfram
% Minimize[{6*x^2 + 6*x*y + 10/3*y^2 , 5/4*x + 3*y <= 3/2, x>=0,y>=0, x+y=1},{x,y}]

%jouet 1
Q = [2 1; 1 2];
e = [4; 5];
Re = 9/2;
sol_exacte = [1; 1]/2;

%jouet 2
% Q = [3 1 0; 1 5 1; 0 1 2];
% e = [1; 2; 3];
% Re = 3;
% sol_exacte = [0; 0; 1];

% ------------ Résolution du problème

param = struct('alpha', @(k) 1, ... #convexite auxiliaire
               'eps', @(k) .1, ... #auxiliaire 
               'beta', 2, ...  #actualisation des prix 
               'kmax', 1000,...
               'seuil', 1e-6);

tic
u = problem_1(Q, e, Re, param);
toc

verif_sol = @(x) struct('J', [x'*Q*x, sol_exacte'*Q*sol_exacte], ...
                        'l', [sum(x) 1], ...
                        'Re',[e'*x, Re], ...
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
hold off;
legend('u_1', 'u_1 exact', 'u_2', 'u_2 exact');
title('Visualisation de la convergence');
ylim([0 1]);