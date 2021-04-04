%% Creation d'un exemple pour le probleme 2 et resolution par PPA

clear 
%----------------------- Données
% commenter et decommenter les exemples selon celui a tester

%exemple 1
Q = [3/2 0; 0 2];
e = [3; 2];
De = 1;

%exemple 2
% Q = [1 1; 1 2];
% e = [1; 2];
% De = 3/2;

%----------------------- Résolution du problème
%Resolution par le probleme auxiliaire : hyperparametres a ajuster selon l'exemple 
param = struct('alpha', @(k) 1.5, ... %convexite auxiliaire
               'eps', @(k) .1, ...    %auxiliaire 
               'beta', 10, ...        %actualisation des prix 
               'kmax', 10000,...
               'seuil', 1e-6);
tic
u = problem2_PPA(Q, e, De, param);
disp(['Temps resolution : ', num2str(toc)])

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

%% Affichage 3D pour les exemples de dimension 2
figure(2);

N = length(e);
if(N~=2)
    disp('Attention : ne traite que les exemples de dimension 2');
    break;
end

%fonction et mesh
fun = @(x,y) [x y]*e;
c = @(x,y) [x,y]*Q*[x;y] <= De & x>=0 & y>=0;

%mesh surface
x=-.5:.002:1;
y=-.5:.002:1;
[X, Y] = meshgrid(x,y);

J = arrayfun(fun, X, Y);
J_not = J;

mask = arrayfun(c, X, Y);
J(~mask) = NaN;
J_not(mask) = NaN;

%mesh ligne
x_line = 0:.01:1;
y_line = 1-x_line;
mask_line = arrayfun(c, x_line, y_line);
J_line = arrayfun(fun, x_line, y_line);
J_line(~mask_line) = NaN;

%affichage surf 
surf(x,y, J', 'EdgeColor', 'None');
hold on;
surf(x,y, J_not', 'EdgeColor','None', 'FaceColor', [.8 .8 .8], 'FaceAlpha', .5); 

%affichage ligne + point
plot3(y_line, x_line, J_line, 'Color', [0    0.7778    0.6111], 'LineWidth', 1.5);
scatter3(u(2, end), u(1, end), fun(u(1, end), u(2, end)), 50, 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0.6 0.2]);
hold off;

%detail affichage
colormap(winter);
title('Affichage de la solution optimale du probleme');
legend('Contraintes inegalité', 'Extérieur domaine', 'Contrainte egalité', 'Solution optimale');
xlabel('u_1'); ylabel('u_2'); zlabel('J(u_1,u_2)');
view([-20 10]);
