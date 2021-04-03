clear 
%----------------------- Données
% commenter et decommenter les exemples selon celui a tester

%exemple 1
% Q = [3/2 0; 0 2];
% e = [3; 2];
% De = 1;

%exemple 2
Q = [1 1; 1 2];
e = [1; 2];
De = 3/2;

%----------------------- Résolution du problème
%Resolution par le probleme auxiliaire
param = struct('alpha', @(k) 1.5, ... #convexite auxiliaire
               'eps', @(k) .1, ... #auxiliaire 
               'beta', 10, ...  #actualisation des prix 
               'kmax', 10000,...
               'seuil', 1e-6);
tic
u = problem2_PPA(Q, e, De, param);
toc

%resolution par le solveur matlab
sol_exacte = problem2_Solver(e, Q, De);

%----------------------- Verification de la solution
str = @(x) num2str(x);
verif_sol = @(x) struct('J', e'*x , ...
                        'J_exa', e'*sol_exacte,...
                        'l', [str(sum(x)) ' == ' str(1)], ...
                        'De',[str(x'*Q*x) ' <= ' str(De)], ...
                        'u', x', ...
                        'u_exa', sol_exacte');

verif_sol(u(:,end))

%% Affichage de la solution
figure(1);

color = ['b' 'r'];
for i = 1:length(color)
    plot(u(i, :), 'Color', color(i));
    hold on;
    plot([0 size(u, 2)], sol_exacte(i)*[1 1], ':','Color', color(i), 'LineWidth', .2);
end
hold off
title('Evolution des u_i en fonction de l`itération');
legend('u_1', 'u_1 exact', 'u_2', 'u_2 exact', 'Location', 'best');
xlabel('itération'); ylabel('u_i');
ylim([0 1]);
