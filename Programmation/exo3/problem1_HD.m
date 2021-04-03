clear all;

N = 10;

Q = diag(4 + (1:N)/2) + diag((N-1):-1:1, 1)/4 + + diag(ones(N-2, 1), 2)/4;
Q = Q + triu(Q, 1)';

e = (1:N)';
Re = ones(1,N)*e/N -1;

min( eig(Q) )

%%
N = 30;

Q = diag(4 + (1:N)/2) + diag((N-1):-1:1, 1)/16 + + diag(ones(N-2, 1), 2)/4;
Q = Q + triu(Q, 1)';

e = (1:N)';
Re = ones(1,N)*e/N +1;

min(eig(Q))
%%

tic
exa = problem1_Solver(Q, e,Re);
toc

param = struct('alpha', @(k) 2, ... #convexite auxiliaire
               'eps', @(k) .1, ... #auxiliaire 
               'beta', 4, ...  #actualisation des prix 
               'kmax', 1000,...
               'seuil', 1e-6);

tic
u = problem1_PPA(Q, e, Re, param);
toc

str = @(x) num2str(x);

test_sol = @(x) struct('J', x'*Q*x, ...
                       'l', [str(sum(x)), ' == ', str(1)], ...
                       'Re',[str(e'*x), ' <= ', str(Re)], ...
                       'u', x');
                     
test_sol(u(:,end))

disp('Comparaison des solutions');
[u(:,end), exa]

%%
figure(2);
col = @(i) [i/N 0 1-i/N];

for i = 1:N
    plot([30 size(u,2)], exa(i)*[1 1], ':', 'Color', col(i)/2, 'LineWidth',1.4);
    hold on;
    plot(u(i, :), 'Color', col(i));
end

legend('Solution PPA', 'Solution point intérieur');

hold off;


