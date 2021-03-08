%Exercice 2:

%% Question 1

%Test des decomposition par les prix, par les quantites et par prediction sur les exemples jouets:

N1=1;

A1=zeros(N1+1,2,2);
A1(1,:,:)=[4 0 ; 0 5];
A1(2,:,:)=[1 0 ; 0 2];

C1=zeros(N1+1,2);
C1(1,:)=[1 2];
C1(2,:)=[1 3];

rho = 0.1;
eps = 10^(-6);
kmax = 4000;

[u1_prix,v1_prix,p1_prix,k1_prix,J1_prix] = DecompositionPrix(N1,A1,C1,rho,eps,kmax);
