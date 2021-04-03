clear
%----------------------- Données
% commenter et decommenter les exemples selon celui a tester

%exemple 1
Q = [2 1; 1 2];
e = [4; 5];
Re = 9/2;
sol_exacte = [1; 1]/2;
%[1.5 .1 10]

%exemple 2
Q = [6 3; 3 10/3];
e = [5/4; 3];
Re = 5/2;
sol_exacte = [1; 9]/10;

%exemple 3
Q = [3 1 0; 1 5 1; 0 1 2];
e = [1; 2; 3];
Re = 1;
sol_exacte = [8; 1; 12]/21;


% ------------ Résolution du problème
%Resolution par le probleme auxiliaire
param = struct('alpha', @(k) 2, ... #convexite auxiliaire
               'eps', @(k) .1, ... #auxiliaire 
               'beta', 4, ...  #actualisation des prix 
               'kmax', 1000,...
               'seuil', 1e-6);

tic
u = problem1_PPA(Q, e, Re, param);
toc

%resolution par le solveur matlab
sol_exacte = problem1_Solver(Q,e,Re);

%----------------------- Verification de la solution
str = @(x) num2str(x);
verif_sol = @(x) struct('J', x'*Q*x , ...
                        'J_exa', sol_exacte'*Q*sol_exacte,...
                        'l', [str(sum(x)) ' == ' str(1)], ...
                        'Re',[str(e'*x) ' >= ' str(Re)], ...
                        'u', x', ...
                        'u_exa', sol_exacte');

verif_sol(u(:,end))

%% Affichage de la solution
figure(1);

N = length(e);
color = winter(N);
for i = 1:N
    plot(u(i, :), 'Color', color(i,:));
    hold on;
    plot([0 size(u, 2)], sol_exacte(i)*[1 1], ':','Color', color(i,:), 'LineWidth', .2);
end
hold off;
title('Evolution des u_i en fonction de l`itération');
legend('u_i', 'u_i exact','Location', 'best');
xlabel('itération'); ylabel('u_i');
ylim([0 1]);