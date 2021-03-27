%Exercice 2:

%% Question 1

%Test des decomposition par les prix, par les quantites et par prediction sur les exemples jouets:

N1=1;

A1=zeros(N1+1,2,2);
A1(1,:,:)=[4 0 ; 0 5];
A1(2,:,:)=[1 0 ; 0 2];

C1=zeros(N1,2);
C1(1,:)=[1 2];

rho = 1;
eps = 10^(-9);
kmax = 100;

[u1_prix,v1_prix,p1_prix,k1_prix,J1_prix] = DecompositionPrix(N1,A1,C1,rho,eps,kmax);
u1_prix
v1_prix

eps = 10^(-9);
rho=0.1;
kmax=100;

[u1_quantites,v1_quantites,w1_quantites,k1_quantites,J1_quantites] = DecompositionQuantites(N1,A1,C1,rho,eps,kmax);
u1_quantites
v1_quantites

eps = 10^(-9);
kmax=100;

[u1_pred_par,v1_pred_par,w1_pred_par,p1_pred_par,k1_pred_par,J1_pred_par] = DecompositionPredictionPar(N1,A1,C1,eps,kmax);
u1_pred_par
v1_pred_par

eps = 10^(-9);
kmax=100;

[u1_pred_seq,v1_pred_seq,w1_pred_seq,p1_pred_seq,k1_pred_seq,J1_pred_seq] = DecompositionPredictionSeq(N1,A1,C1,eps,kmax);
u1_pred_seq
v1_pred_seq