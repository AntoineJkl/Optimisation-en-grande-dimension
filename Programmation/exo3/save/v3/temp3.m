clear all;

N = 10;

Q = diag(4 + (1:N)/2) + diag((N-1):-1:1, 1)/4 + + diag(ones(N-2, 1), 2)/4;
Q = Q + triu(Q, 1)';

e = (1:N)';
De = ones(1,N)*Q*ones(N,1)/N^2 +1;

min( eig(Q) )

%%
N = 30;

Q = diag(4 + (1:N)/2) + diag((N-1):-1:1, 1)/16 + + diag(ones(N-2, 1), 2)/4;
Q = Q + triu(Q, 1)';

e = (1:N)';
De = ones(1,N)*Q*ones(N,1)/N^2 +1;

min(eig(Q))
%%
x0 = ones(N,1)/N;

tic
exa = ResolutionExactProbleme2(N,x0,e, Q, De);
toc

param = struct('alpha', @(k) 4, ... #convexite auxiliaire
               'eps', @(k) .1, ... #auxiliaire 
               'beta', 4, ...  #actualisation des prix 
               'kmax', 1000,...
               'seuil', 1e-6);

tic
u = problem_2(Q, e, De, param);
toc

[u(:,end), exa]


str = @(x) num2str(x);

test_sol = @(x) struct('J', e'*x, ...
                       'l', [str(sum(x)), ' == ', str(1)], ...
                       'Re',[str(x'*Q*x), ' <= ', str(De)], ...
                       'u', x');

disp('---------------------------------------');
disp('Comparaison des solutions');
test_sol(u(:,end))


%%
figure(2);
col = @(i) [i/N 0 1-i/N];

for i = 1:N
    plot([0 size(u,2)], exa(i)*[1 1], ':', 'Color', col(i)/2, 'LineWidth',1.4);
    hold on;
    plot(u(i, :), 'Color', col(i));
end

ylim([-.01 .4]);

legend('Solution point intérieur', 'Solution PPA');

hold off;
