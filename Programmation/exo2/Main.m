%Exercice 2:

%% Question 1

%% Test des decomposition par les prix, par les quantites et par prediction sur les exemples jouets ( programmation Uzawa ou Arrow ) :
addpath('.\sequentielle');

N1=1;

A1=zeros(N1+1,2,2);
A1(1,:,:)=[4 0 ; 0 5];
A1(2,:,:)=[1 0 ; 0 2];

C1=zeros(N1,2);
C1(1,:)=[1 2];

param_prix = struct('rho', 1, ...
                    'eps', 10^(-3), ...
                    'algo','arrow',...
                    'kmax', 100);
disp('Solution par la décomposition par les prix');
[u1_prix,v1_prix,~,k1_prix,~,u_stock_prix] = DecompositionPrix(N1,A1,C1,param_prix);
u1_prix
v1_prix



param_quantite1 = struct('rho', 0.4, ...
                    'eps', 10^(-5), ...
                    'algo','uzawa',...
                    'kmax', 100);

disp('Solution par la décomposition par les quantités');
[u1_quantites,v1_quantites,~,k1_quantites,J,u_stock_quantites1] = DecompositionQuantites(N1,A1,C1,param_quantite1);
u1_quantites
v1_quantites


param_quantite2 = struct('rho', 0.3, ...
                    'algo','uzawa',...
                    'eps', 10^(-5), ...
                    'kmax', 100);
                
disp('Solution par la décomposition par les quantités avec modifications');
[u1_quantites2,v1_quantites2,~,k1_quantites2,~,u_stock_quantites2] = DecompositionQuantites2(N1,A1,C1,param_quantite2);
u1_quantites2
v1_quantites2

param_pred_par = struct('eps', 10^(-5), ...
                    'algo','arrow',...
                    'beta', 0.5, ...
                    'gamma', 0.5 ,...
                    'kmax', 100);
                
disp('Solution par la décomposition par prédiction parallèle');
[u1_pred_par,v1_pred_par,~,~,k1_pred_par,~,u_stock_pred_par] = DecompositionPredictionPar(N1,A1,C1,param_pred_par);
u1_pred_par
v1_pred_par

param_pred_seq = struct('eps', 10^(-5), ...
                    'algo','arrow',...
                    'beta', 0.5, ...
                    'gamma', 1 ,...
                    'kmax', 100);

disp('Solution par la décomposition par prédiction séquentielle 1');
[u1_pred_seq,v1_pred_seq,~,~,k1_pred_seq,~,u_stock_pred_seq] = DecompositionPredictionSeq(N1,A1,C1,param_pred_seq);
u1_pred_seq
v1_pred_seq


param_pred_seq2 = struct('eps', 10^(-5), ...
                    'algo','arrow',...
                    'beta', 0.5, ...
                    'gamma', 1 ,...
                    'kmax', 100);

disp('Solution par la décomposition par prédiction séquentielle 2');
[u1_pred_seq2,v1_pred_seq2,~,~,k1_pred_seq2,~,u_stock_pred_seq2] = DecompositionPredictionSeq2(N1,A1,C1,param_pred_seq2);
u1_pred_seq2
v1_pred_seq2


%% Test des decomposition par les quantites sur les exemples jouets ( programmation explicite ) :
addpath('.\sequentielle');

param_quantite1 = struct('rho', 0.4, ...
                    'eps', 10^(-5), ...
                    'algo','explicite',...
                    'kmax', 100);

disp('Solution par la décomposition par les quantités');
[u1_quantites,v1_quantites,~,k1_quantites,~,u_stock_quantites1] = DecompositionQuantites(N1,A1,C1,param_quantite1);
u1_quantites
v1_quantites


param_quantite2 = struct('rho', 0.3, ...
                    'algo','explicite',...
                    'eps', 10^(-5), ...
                    'kmax', 100);
                
disp('Solution par la décomposition par les quantités avec modifications');
[u1_quantites2,v1_quantites2,~,k1_quantites2,~,u_stock_quantites2] = DecompositionQuantites2(N1,A1,C1,param_quantite2);
u1_quantites2
v1_quantites2

%% Test des decomposition par les prix, par les quantites et par prediction sur les exemples jouets ( programmation méthode point intérieur ) :
addpath('.\sequentielle');

param_prix = struct('rho', 1, ...
                    'eps', 10^(-3), ...
                    'kmax', 100);
disp('Solution par la décomposition par les prix');
[u1_prix,v1_prix,~,k1_prix,~,u_stock_prix] = DecompositionPrix(N1,A1,C1,param_prix);
u1_prix
v1_prix



param_quantite1 = struct('rho', 0.4, ...
                    'eps', 10^(-5), ...
                    'kmax', 100);

disp('Solution par la décomposition par les quantités');
[u1_quantites,v1_quantites,w1_quantites,k1_quantites,~,u_stock_quantites1] = DecompositionQuantites(N1,A1,C1,param_quantite1);
u1_quantites
v1_quantites


param_quantite2 = struct('rho', 0.3, ...
                    'eps', 10^(-5), ...
                    'kmax', 100);
                
disp('Solution par la décomposition par les quantités avec modifications');
[u1_quantites2,v1_quantites2,~,k1_quantites2,~,u_stock_quantites2] = DecompositionQuantites2(N1,A1,C1,param_quantite2);
u1_quantites2
v1_quantites2

param_pred_par = struct('eps', 10^(-5), ...
                    'beta', 0.5, ...
                    'gamma', 0.5 ,...
                    'kmax', 100);
                
disp('Solution par la décomposition par prédiction parallèle');
[u1_pred_par,v1_pred_par,~,~,k1_pred_par,~,u_stock_pred_par] = DecompositionPredictionPar(N1,A1,C1,param_pred_par);
u1_pred_par
v1_pred_par

param_pred_seq = struct('eps', 10^(-5), ...
                    'beta', 0.5, ...
                    'gamma', 1 ,...
                    'kmax', 100);

disp('Solution par la décomposition par prédiction séquentielle 1');
[u1_pred_seq,v1_pred_seq,~,~,k1_pred_seq,~,u_stock_pred_seq] = DecompositionPredictionSeq(N1,A1,C1,param_pred_seq);
u1_pred_seq
v1_pred_seq


param_pred_seq2 = struct('eps', 10^(-5), ...
                    'beta', 0.5, ...
                    'gamma', 1 ,...
                    'kmax', 100);

disp('Solution par la décomposition par prédiction séquentielle 2');
[u1_pred_seq2,v1_pred_seq2,~,~,k1_pred_seq2,~,u_stock_pred_seq2] = DecompositionPredictionSeq2(N1,A1,C1,param_pred_seq2);
u1_pred_seq2
v1_pred_seq2


%% Affichage question 1 prix

u11=u_stock_prix(1:k1_prix,1,1); u21=u_stock_prix(1:k1_prix,2,1);
u12=u_stock_prix(1:k1_prix,1,2); u22=u_stock_prix(1:k1_prix,2,2);

figure
hold on;
plot(1:k1_prix,u11,1:k1_prix,u12);
plot([0,k1_prix],[1 1],'--ro');
plot([0, k1_prix],[0.5263 0.5263],'--ko');
legend('u11 prix','u12 prix','u11 theo','u12 theo');
xlabel('Iterations');
title('Convergence algorithme de décomposition par les prix');
hold off;

figure
plot(1:k1_prix,abs(u21-u11),1:k1_prix,abs(u22-u12));
legend('prix dim 1','prix dim 2');
xlabel('Iterations');
title('Erreur d egalite au niveau de la contrainte couplante (prix)')

%% Affichage question 1 quantites

u111=u_stock_quantites1(1:k1_quantites,1,1); u211=u_stock_quantites1(1:k1_quantites,2,1);
u121=u_stock_quantites1(1:k1_quantites,1,2); u221=u_stock_quantites1(1:k1_quantites,2,2);

u112=u_stock_quantites2(1:k1_quantites2,1,1); u212=u_stock_quantites2(1:k1_quantites2,2,1);
u122=u_stock_quantites2(1:k1_quantites2,1,2); u222=u_stock_quantites2(1:k1_quantites2,2,2);

figure
hold on;
plot(1:k1_quantites,u111,1:k1_quantites,u121);
plot(1:k1_quantites2,u112,1:k1_quantites2,u122);
plot([0,max(k1_quantites,k1_quantites2)],[1 1],'--ro');
plot([0,max(k1_quantites,k1_quantites2)],[0.5263 0.5263],'--ko');
plot([k1_quantites,k1_quantites],[0 2],'--r','Color',[0.5, 0.3, 0.15]);
plot([k1_quantites2,k1_quantites2],[0 2],'--g');
xlabel('Iterations');
legend('u11 quantites','u12 quantites','u11 quantites2','u12 quantites1','u11 theo','u12 theo');
title('Convergence algorithme de décomposition par les quantites');
hold off;

figure
hold on;
plot(1:k1_quantites,abs(u211-u111),1:k1_quantites,abs(u221-u121));
plot(1:k1_quantites2,abs(u212-u112),1:k1_quantites2,abs(u222-u122));
ylim([0,1]);
xlabel('Iterations');
legend('quantites dim 1','quantites dim 2','quantites2 dim 1','quantites2 dim 2');
title('Erreur d egalite au niveau de la contrainte couplante (quantite)')
hold off;

%% Affichage question 1 prediction

u111=u_stock_pred_par(1:k1_pred_par,1,1); u211=u_stock_pred_par(1:k1_pred_par,2,1);
u121=u_stock_pred_par(1:k1_pred_par,1,2); u221=u_stock_pred_par(1:k1_pred_par,2,2);

u112=u_stock_pred_seq(1:k1_pred_seq,1,1); u212=u_stock_pred_seq(1:k1_pred_seq,2,1);
u122=u_stock_pred_seq(1:k1_pred_seq,1,2); u222=u_stock_pred_seq(1:k1_pred_seq,2,2);

u113=u_stock_pred_seq2(1:k1_pred_seq2,1,1); u213=u_stock_pred_seq2(1:k1_pred_seq2,2,1);
u123=u_stock_pred_seq2(1:k1_pred_seq2,1,2); u223=u_stock_pred_seq2(1:k1_pred_seq2,2,2);

figure
hold on;
plot(1:k1_pred_par,u111,1:k1_pred_par,u121);
plot(1:k1_pred_seq,u112,1:k1_pred_seq,u122);
plot(1:k1_pred_seq2,u113,1:k1_pred_seq2,u123);
plot([0,max(k1_pred_seq,k1_pred_par)],[1 1],'--ro');
plot([0,max(k1_pred_seq,k1_pred_par)],[0.5263 0.5263],'--ko');
plot([k1_pred_par,k1_pred_par],[0 2],'--r','Color',[0.5, 0.3, 0.15]);
plot([k1_pred_seq,k1_pred_seq],[0 2],'--g');
xlabel('Iterations');
legend('u11 pred par','u12 pred par','u11 pred seq','u12 pred seq','u11 pred seq2','u12 pred seq2','u11 theo','u12 theo');
title('Convergence algorithme de décomposition par prédiction');
hold off;

figure
hold on;
plot(1:k1_pred_par,abs(u211-u111),1:k1_pred_par,abs(u221-u121));
plot(1:k1_pred_seq,abs(u212-u112),1:k1_pred_seq,abs(u222-u122));
plot(1:k1_pred_seq,abs(u213-u113),1:k1_pred_seq,abs(u223-u123));
ylim([0,1]);
xlabel('Iterations');
legend('pred par dim 1','pred par dim 2','pred seq dim 1','pred seq dim 2','pred seq dim 1','pred seq dim 2');
title('Erreur d egalite au niveau de la contrainte couplante (prédiction)')
hold off;

%% Question 2 : Différence 1 - Relaxation ou Non Relaxation
N1=1;

A1=zeros(N1+1,2,2);
A1(1,:,:)=[4 0 ; 0 5];
A1(2,:,:)=[1 0 ; 0 2];

C1=zeros(N1,2);
C1(1,:)=[1 2];



addpath('.\sequentielle');

param_pred_par = struct('eps', 10^(-2), ...
                    'beta', 1, ...
                    'gamma', 1 ,...
                    'kmax', 30);

disp('Solution par la décomposition par prédiction parallèle');
[u2_pred_par,v2_pred_par,w2_pred_par,p2_pred_par,k2_pred_par,J2_pred_par,u_stock_pred_par] = DecompositionPredictionPar(N1,A1,C1,param_pred_par);
u2_pred_par

param_pred_seq = struct('eps', 10^(-2), ...
                    'beta', 1, ...
                    'gamma', 1 ,...
                    'kmax', 30);

disp('Solution par la décomposition par prédiction séquentielle');
[u2_pred_seq,v2_pred_seq,w2_pred_seq,p2_pred_seq,k2_pred_seq,J2_pred_seq,u_stock_pred_seq] = DecompositionPredictionSeq(N1,A1,C1,param_pred_seq);
u2_pred_seq

param_pred_seq2 = struct('eps', 10^(-2), ...
                    'beta', 1, ...
                    'gamma', 1 ,...
                    'kmax', 30);

disp('Solution par la décomposition par prédiction séquentielle 2');
[u2_pred_seq2,v2_pred_seq2,w2_pred_seq2,p2_pred_seq2,k2_pred_seq2,J2_pred_seq2,u_stock_pred_seq2] = DecompositionPredictionSeq2(N1,A1,C1,param_pred_seq2);
u2_pred_seq2

%% Affichage Différence 1 - Relaxation ou Non Relaxation

u111=u_stock_pred_par(1:k2_pred_par,1,1); u211=u_stock_pred_par(1:k2_pred_par,2,1);
u121=u_stock_pred_par(1:k2_pred_par,1,2); u221=u_stock_pred_par(1:k2_pred_par,2,2);

u112=u_stock_pred_seq(1:k2_pred_seq,1,1); u212=u_stock_pred_seq(1:k2_pred_seq,2,1);
u122=u_stock_pred_seq(1:k2_pred_seq,1,2); u222=u_stock_pred_seq(1:k2_pred_seq,2,2);

u113=u_stock_pred_seq2(1:k2_pred_seq2,1,1); u213=u_stock_pred_seq2(1:k2_pred_seq2,2,1);
u123=u_stock_pred_seq2(1:k2_pred_seq2,1,2); u223=u_stock_pred_seq2(1:k2_pred_seq2,2,2);

figure
hold on;
plot(1:k2_pred_par,u111,1:k2_pred_par,u121);
plot([0,max(k2_pred_seq,k2_pred_par)],[1 1],'--ro');
plot([0,max(k2_pred_seq,k2_pred_par)],[0.5263 0.5263],'--ko');
legend('u11 pred par','u12 pred par','u11 theo','u12 theo');
title('Convergence algorithme de décomposition par prédiction parralèle');
hold off;

figure
hold on;
plot(1:k2_pred_seq,u112,1:k2_pred_seq,u122);
plot([0,k2_pred_seq],[1 1],'--ro');
plot([0,k2_pred_seq],[0.5263 0.5263],'--ko');
legend('u11 pred seq','u12 pred seq','u11 theo','u12 theo');
title('Convergence algorithme de décomposition par prédiction séquentielle');
hold off;

figure
hold on;
plot(1:k2_pred_seq2,u113,1:k2_pred_seq2,u123);
plot([0,k2_pred_seq2],[1 1],'--ro');
plot([0,k2_pred_seq2],[0.5263 0.5263],'--ko');
legend('u11 pred seq2','u12 pred seq2','u11 theo','u12 theo');
title('Convergence algorithme de décomposition par prédiction séquentielle 2');
hold off;

%% Question 2 : Différence 2 - Relaxation sur les prix

N1=1;

A1=zeros(N1+1,2,2);
A1(1,:,:)=[4 0 ; 0 5];
A1(2,:,:)=[1 0 ; 0 2];

C1=zeros(N1,2);
C1(1,:)=[1 2];



addpath('.\sequentielle');

param_pred_par = struct('eps', 10^(-2), ...
                    'beta', 0.5, ...
                    'gamma', 1 ,...
                    'kmax', 30);

disp('Solution par la décomposition par prédiction parallèle');
[u2_pred_par,v2_pred_par,w2_pred_par,p2_pred_par,k2_pred_par,J2_pred_par,u_stock_pred_par] = DecompositionPredictionPar(N1,A1,C1,param_pred_par);
u2_pred_par

param_pred_seq = struct('eps', 10^(-2), ...
                    'beta', 0.5, ...
                    'gamma', 1 ,...
                    'kmax', 30);

disp('Solution par la décomposition par prédiction séquentielle');
[u2_pred_seq,v2_pred_seq,w2_pred_seq,p2_pred_seq,k2_pred_seq,J2_pred_seq,u_stock_pred_seq] = DecompositionPredictionSeq(N1,A1,C1,param_pred_seq);
u2_pred_seq

param_pred_seq2 = struct('eps', 10^(-2), ...
                    'beta', 0.5, ...
                    'gamma', 1 ,...
                    'kmax', 30);

disp('Solution par la décomposition par prédiction séquentielle 2');
[u2_pred_seq2,v2_pred_seq2,w2_pred_seq2,p2_pred_seq2,k2_pred_seq2,J2_pred_seq2,u_stock_pred_seq2] = DecompositionPredictionSeq2(N1,A1,C1,param_pred_seq2);
u2_pred_seq2

%% Affichage Différence 2 - Relaxation sur les prix

u111=u_stock_pred_par(1:k2_pred_par,1,1); u211=u_stock_pred_par(1:k2_pred_par,2,1);
u121=u_stock_pred_par(1:k2_pred_par,1,2); u221=u_stock_pred_par(1:k2_pred_par,2,2);

u112=u_stock_pred_seq(1:k2_pred_seq,1,1); u212=u_stock_pred_seq(1:k2_pred_seq,2,1);
u122=u_stock_pred_seq(1:k2_pred_seq,1,2); u222=u_stock_pred_seq(1:k2_pred_seq,2,2);

u113=u_stock_pred_seq2(1:k2_pred_seq2,1,1); u213=u_stock_pred_seq2(1:k2_pred_seq2,2,1);
u123=u_stock_pred_seq2(1:k2_pred_seq2,1,2); u223=u_stock_pred_seq2(1:k2_pred_seq2,2,2);

figure
hold on;
plot(1:k2_pred_par,u111,1:k2_pred_par,u121);
plot([0,max(k2_pred_seq,k2_pred_par)],[1 1],'--ro');
plot([0,max(k2_pred_seq,k2_pred_par)],[0.5263 0.5263],'--ko');
legend('u11 pred par','u12 pred par','u11 theo','u12 theo');
title('Convergence algorithme de décomposition par prédiction parralèle');
hold off;

figure
hold on;
plot(1:k2_pred_seq,u112,1:k2_pred_seq,u122);
plot([0,k2_pred_seq],[1 1],'--ro');
plot([0,k2_pred_seq],[0.5263 0.5263],'--ko');
legend('u11 pred seq','u12 pred seq','u11 theo','u12 theo');
title('Convergence algorithme de décomposition par prédiction séquentielle');
hold off;

figure
hold on;
plot(1:k2_pred_seq2,u113,1:k2_pred_seq2,u123);
plot([0,k2_pred_seq2],[1 1],'--ro');
plot([0,k2_pred_seq2],[0.5263 0.5263],'--ko');
legend('u11 pred seq2','u12 pred seq2','u11 theo','u12 theo');
title('Convergence algorithme de décomposition par prédiction séquentielle 2');
hold off;


%% Question 2 : Différence 3 - Non respect de la contrainte couplante avec la relaxation sur les allocations

N1=1;

A1=zeros(N1+1,2,2);
A1(1,:,:)=[4 0 ; 0 5];
A1(2,:,:)=[1 0 ; 0 2];

C1=zeros(N1,2);
C1(1,:)=[1 2];

addpath('.\sequentielle');

param_pred_par = struct('eps', 10^(-3), ...
                    'beta', 0.5, ...
                    'gamma', 0.5 ,...
                    'kmax', 30);

disp('Solution par la décomposition par prédiction parallèle');
[u2_pred_par,v2_pred_par,w2_pred_par,p2_pred_par,k2_pred_par,J2_pred_par,u_stock_pred_par] = DecompositionPredictionPar(N1,A1,C1,param_pred_par);
u2_pred_par

param_pred_seq = struct('eps', 10^(-3), ...
                    'beta', 0.5, ...
                    'gamma', 1 ,...
                    'kmax', 30);

disp('Solution par la décomposition par prédiction séquentielle');
[u2_pred_seq,v2_pred_seq,w2_pred_seq,p2_pred_seq,k2_pred_seq,J2_pred_seq,u_stock_pred_seq] = DecompositionPredictionSeq(N1,A1,C1,param_pred_seq);
u2_pred_seq



%% Affichage Différence 3 - Non respect de la contrainte couplante avec la relaxation sur les allocations

u111=u_stock_pred_par(1:k2_pred_par,1,1); u211=u_stock_pred_par(1:k2_pred_par,2,1);
u121=u_stock_pred_par(1:k2_pred_par,1,2); u221=u_stock_pred_par(1:k2_pred_par,2,2);

u112=u_stock_pred_seq(1:k2_pred_seq,1,1); u212=u_stock_pred_seq(1:k2_pred_seq,2,1);
u122=u_stock_pred_seq(1:k2_pred_seq,1,2); u222=u_stock_pred_seq(1:k2_pred_seq,2,2);

u113=u_stock_pred_seq2(1:k2_pred_seq2,1,1); u213=u_stock_pred_seq2(1:k2_pred_seq2,2,1);
u123=u_stock_pred_seq2(1:k2_pred_seq2,1,2); u223=u_stock_pred_seq2(1:k2_pred_seq2,2,2);

figure
hold on;
plot(1:k2_pred_par,u111,1:k2_pred_par,u121);
plot(1:k2_pred_seq,u112,1:k2_pred_seq,u122);
plot([0,max(k2_pred_seq,k2_pred_par)],[1 1],'--ro');
plot([0,max(k2_pred_seq,k2_pred_par)],[0.5263 0.5263],'--ko');
plot([k2_pred_par,k2_pred_par],[0 2],'--r','Color',[0.5, 0.3, 0.15]);
plot([k2_pred_seq,k2_pred_seq],[0 2],'--g');
legend('u11 pred par','u12 pred par','u11 pred seq','u12 pred seq','u11 theo','u12 theo');
title('Convergence algorithme de décomposition par prédiction parralèle');
hold off;

figure
hold on;
plot(1:k2_pred_par,abs(u211-u111),1:k2_pred_par,abs(u221-u121));
plot(1:k2_pred_seq,abs(u212-u112),1:k2_pred_seq,abs(u222-u122));
ylim([0,1]);
legend('pred par dim 1','pred par dim 2','pred seq dim 1','pred seq dim 2','pred seq dim 1','pred seq dim 2');
title('Erreur d egalite au niveau de la contrainte couplante (prédiction)')
hold off;

%% Question 2 : Différence 4 - Temps d'execution et nombre d'itérations en fonction de N 

N_tab = 1:1:5;
k_tab_par = []; k_tab_seq = []; k_tab_seq2 = []; k_tab_seq22 = [] ;
t_tab_par = []; t_tab_seq = []; t_tab_seq2 = []; t_tab_seq22 = [];
J_tab_par = []; J_tab_seq = []; J_tab_seq2 = []; J_tab_seq22 = [];

disp('Etude temps d execution et nombre d iterations en fonction de N');
for N = N_tab
    disp(['N = ' num2str(N)]);
    [A,C] = CreationInstance(N,1,3,1,3);

    addpath('.\sequentielle');
    
    
    %Decomposition par prédiction parallèle
    param_pred_par = struct('eps', 10^(-1), ...
                    'beta', 0.5, ...
                    'gamma', 0.5 ,...
                    'kmax', 100);
    tic
    [~,~,~,~,k1_pred_par,J_par] = DecompositionPredictionPar(N,A,C,param_pred_par);
    t_tab_par=[t_tab_par,toc];
    k_tab_par=[k_tab_par , k1_pred_par];
    J_tab_par=[J_tab_par,J_par];
    
    %Decomposition par prédiction séquentielle
    param_pred_seq = struct('eps', 10^(-1), ...
                    'beta', 0.5, ...
                    'gamma', 0.5 ,...
                    'kmax', 100);
    tic
    [~,~,~,~,k1_pred_seq,J_seq] = DecompositionPredictionSeq(N,A,C,param_pred_seq);
    t_tab_seq=[t_tab_seq,toc];
    k_tab_seq=[k_tab_seq , k1_pred_seq];
    J_tab_seq=[J_tab_seq,J_seq];
    
    %Decomposition par prédiction séquentielle
    param_pred_seq = struct('eps', 10^(-1), ...
                    'beta', 0.5, ...
                    'gamma', 1 ,...
                    'kmax', 100);
    tic
    [~,~,~,~,k1_pred_seq22,J_seq22] = DecompositionPredictionSeq(N,A,C,param_pred_seq);
    t_tab_seq22=[t_tab_seq22,toc];
    k_tab_seq22=[k_tab_seq22 , k1_pred_seq22];
    J_tab_seq22=[J_tab_seq22,J_seq22];
    
    
    %Decomposition par prédiction séquentielle 2
    param_pred_seq2 = struct('eps', 10^(-1), ...
                    'beta', 0.5, ...
                    'gamma', 0.5 ,...
                    'kmax', 100);
                
    tic
    [~,~,~,~,k1_pred_seq2,J_seq2] = DecompositionPredictionSeq2(N,A,C,param_pred_seq2);
    t_tab_seq2=[t_tab_seq2,toc];
    k_tab_seq2=[k_tab_seq2 , k1_pred_seq2];
    J_tab_seq2=[J_tab_seq2,J_seq2];

end

%% Affichage : Différence 4 - Temps d'execution et nombre d'itérations en fonction de N 

figure;
plot(N_tab,t_tab_par,N_tab,t_tab_seq,N_tab,t_tab_seq22,N_tab,t_tab_seq2);
legend('parallele','sequentielle 1 gamma=0.5','sequentielle 1 gamma=1','sequentielle 2');
title({'Temps d execution en fonction de N des'  'différentes décompositions par prédiction'});
xlabel('N'); ylabel('Temps en s');

figure;
hold on;
plot(N_tab,k_tab_par,N_tab,k_tab_seq,N_tab,k_tab_seq22,N_tab,k_tab_seq2);
plot(N_tab,[101,101,101,101,101],'--ro');
legend('parallele','sequentielle 1 gamma=0.5','sequentielle 1 gamma=1','sequentielle 2');
title({'Nombre d iteration en fonction de N des'  'différentes décompositions par prédiction'});
xlabel('N'); ylabel('Iterations');
hold off;

figure;
plot(N_tab,J_tab_par,N_tab,J_tab_seq,N_tab,J_tab_seq22,N_tab,J_tab_seq2);
legend('parallele','sequentielle 1 gamma=0.5','sequentielle 1 gamma=1','sequentielle 2');
title({'Valeur objective finale en fonction de N des'  'différentes décompositions par prédiction'});
xlabel('N'); ylabel('Valeur de la fonction objectrive');


%% Question 3

%% 3.1 - Etude du nombre d'iteration et du temps d'execution en fonction de N
N_tab = 20:20:200;
k_tab_prix = []; k_tab_quantite = [];
t_tab_prix = []; t_tab_quantite = [];
J_tab_prix = [] ; J_tab_quantite = [];
disp('Etude temps d execution et nombre d iterations en fonction de N');
for N = N_tab
    disp(['N = ' num2str(N)]);
    
    %Creation de l'instance
    addpath('.');
    [A,C] = CreationInstance(N,1,2,0,1);
    
    addpath('.\parallele');
    %Decomposition par les prix
    param_prix = struct('rho', 0.1, ...
                    'eps', 10^(-2), ...
                    'kmax', 1000);
    tic
    [~,~,~,k1_prix,J1_prix] = DecompositionPrix_par(N,A,C,param_prix);
    t_tab_prix = [t_tab_prix,toc];
    k_tab_prix = [k_tab_prix , k1_prix];
    J_tab_prix = [J_tab_prix,J1_prix];
    
    %Décomposition par les quantites
    param_quantite2 = struct('rho', 0.1, ...
                        'algo' , 'explicite' , ...
                        'eps', 10^(-2), ...
                        'kmax', 1000);
    tic
    [~,~,~,k1_quantites2,J1_quantites2] = DecompositionQuantites2_par(N,A,C,param_quantite2);
    t_tab_quantite = [t_tab_quantite,toc];
    k_tab_quantite = [k_tab_quantite , k1_quantites2];
    J_tab_quantite = [J_tab_quantite,J1_quantites2];

end

%% Affichage 3.1
figure;
plot(N_tab,t_tab_prix,N_tab,t_tab_quantite);
legend('prix','quantites2');
title('Temps d execution en fonction de N');
xlabel('N'); ylabel('Temps en s');

figure;
plot(N_tab,k_tab_prix,N_tab,k_tab_quantite);
legend('prix','quantites2');
title('Nombre d iterations en fonction de N');
xlabel('N'); ylabel('Iterations');

figure;
plot(N_tab,J_tab_prix,N_tab,J_tab_quantite);
legend('prix','quantites2');
title('Valeur objective en fonction de N');
xlabel('N'); ylabel('Valeur objective');

%% 3.2 - Etude du nombre d'iteration et du temps d'execution en fonction de la précision

%Creation de l'instance
N=10;
[A,C] = CreationInstance(N,1,2,0,1);

k_tab_prix = []; k_tab_quantite = [];
t_tab_prix = []; t_tab_quantite = [];
J_tab_prix = [] ; J_tab_quantite = [];

pho_tab = [0.00001 , 0.00003 , 0.0001 , 0.0003 , 0.001 , 0.003 , 0.01 , 0.03 , 0.1 , 0.3 , 1 , 3];

disp('Etude temps d execution et nombre d iterations en fonction de N');
for pho = pho_tab
    disp(['pho = ' num2str(pho)]);
    
    addpath('.\parallele');
    %Decomposition par les prix
    param_prix = struct('rho', pho, ...
                    'eps', 10^(-2), ...
                    'kmax', 1000);
    tic
    [~,~,~,k1_prix,J1_prix] = DecompositionPrix_par(N,A,C,param_prix);
    t_tab_prix = [t_tab_prix,toc];
    k_tab_prix = [k_tab_prix , k1_prix];
    J_tab_prix = [J_tab_prix,J1_prix];
    
    %Décomposition par les quantites
    param_quantite2 = struct('rho', pho, ...
                        'algo' , 'explicite' , ...
                        'eps', 10^(-2), ...
                        'kmax', 1000);
    tic
    [~,~,~,k1_quantites2,J1_quantites2] = DecompositionQuantites2_par(N,A,C,param_quantite2);
    t_tab_quantite = [t_tab_quantite,toc];
    k_tab_quantite = [k_tab_quantite , k1_quantites2];
    J_tab_quantite = [J_tab_quantite,J1_quantites2];

end

%% Affichage 3.2

figure;
plot(pho_tab,t_tab_prix,pho_tab,t_tab_quantite);
legend('prix','quantites2');
title('Temps d execution en fonction de pho');
xlabel('pho'); ylabel('Temps en s');

figure;
plot(pho_tab,k_tab_prix,pho_tab,k_tab_quantite);
legend('prix','quantites2');
title('Nombre d iterations en fonction de pho');
xlabel('pho'); ylabel('Iterations');

figure;
plot(pho_tab,J_tab_prix,pho_tab,J_tab_quantite);
legend('prix','quantites2');
title('Valeur objective en fonction de pho');
xlabel('pho'); ylabel('Valeur objective');

%% 3.3 - Etude du nombre d'iteration et du temps d'execution en fonction du pas de la descente

%Creation de l'instance
N=10;
[A,C] = CreationInstance(N,1,2,0,1);

k_tab_prix = []; k_tab_quantite = [];
t_tab_prix = []; t_tab_quantite = [];
J_tab_prix = [] ; J_tab_quantite = [];
residus_prix=[] ; residus_quantite = [];

eps_tab = [1 , 10^(-1) , 10^(-2) , 10^(-3) , 10^(-4) , 10^(-5) , 10^(-6) , 10^(-7) , 10^(-8) , 10^(-9) , 10^(-10)];

disp('Etude temps d execution et nombre d iterations en fonction de N');
for eps = eps_tab
    disp(['eps = ' num2str(eps)]);
    
    addpath('.\parallele');
    %Decomposition par les prix
    param_prix = struct('rho', 0.1, ...
                    'eps', eps, ...
                    'kmax', 1000);
    tic
    [u_p,~,~,k1_prix,J1_prix] = DecompositionPrix_par(N,A,C,param_prix);
    residus_prix = [residus_prix, (sum(u_p(1:N,:),1)-u_p(N+1,:))' ];
    t_tab_prix = [t_tab_prix,toc];
    k_tab_prix = [k_tab_prix , k1_prix];
    J_tab_prix = [J_tab_prix,J1_prix];
    
    %Décomposition par les quantites
    param_quantite2 = struct('rho', 0.1, ...
                        'algo' , 'explicite' , ...
                        'eps', eps, ...
                        'kmax', 1000);
    tic
    [u_q,~,~,k1_quantites2,J1_quantites2] = DecompositionQuantites2_par(N,A,C,param_quantite2);
    residus_quantite = [residus_quantite, (sum(u_q(1:N,:),1)-u_q(N+1,:))' ];
    t_tab_quantite = [t_tab_quantite,toc];
    k_tab_quantite = [k_tab_quantite , k1_quantites2];
    J_tab_quantite = [J_tab_quantite,J1_quantites2];
    
end


%% Affichage 3.3

figure;
plot(eps_tab,t_tab_prix,eps_tab,t_tab_quantite);
legend('prix','quantites2');
title('Temps d execution en fonction de la précision');
xlabel('eps'); ylabel('Temps en s');

figure;
plot(eps_tab,k_tab_prix,eps_tab,k_tab_quantite);
legend('prix','quantites2');
title('Nombre d iterations en fonction de la précision');
xlabel('eps'); ylabel('Iterations');

figure;
plot(eps_tab,J_tab_prix,eps_tab,J_tab_quantite);
legend('prix','quantites2');
title('Valeur objective en fonction de la précision');
xlabel('eps'); ylabel('Valeur objective');

figure;
plot(eps_tab,residus_prix(1,:),eps_tab,residus_prix(2,:),eps_tab,residus_quantite(1,:),eps_tab,residus_quantite(2,:));
legend('prix dim1','prix dim2','quantites2 dim1','quantites dim2');
title({'Erreur au niveau de la contrainte' , 'couplante en fonction de la précision'});
xlabel('eps'); ylabel('Erreur');
