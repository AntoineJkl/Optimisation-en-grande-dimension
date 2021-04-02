%Exercice 1:

%% Question 1

%Representation graphique de l'ensemble des solutions admissibles:

%Affichage des droites delimitant le domaine
f1 = @(x) - 1/2*x; f2 = @(x) x*0;
x1 = -2:0.01:2; y1 = f1(x1);
x2 = -2:0.01:2; y2 = f2(x2); 
plot(x1,y1,x2,y2);
hold on;

%Coloration de l'ensemble admissible
ind = @(u1,u2) 1*((u1 + 2*u2 <= 0) & (u2 <= 0));
x_value = -2:.01:2; y_value = -2:.01:2;
[X,Y] = meshgrid(x_value,y_value);
z = ind(X,Y);
objet = pcolor(x_value,y_value,z);
set(objet,'EdgeColor','none');
set(objet,'facealpha',0.3);
colormap(winter);
hold off

title({'Ensemble des solutions','admissibles de J(u)'});
xlabel('u1');
ylabel('u2');
legend('u1 + 2u2 = 0','u2 = 0');






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 2

%Test des decomposition par les prix, par les quantites et par prediction sur les exemples jouets:

%Pour N = 2: 
%   Solution: u* = (1;-1) et J(u*) = -1
N1 = 2;
[A1,b1,C1] = CreateInstance(N1);
rho = 0.1;
eps = 10^(-6);
kmax = 10000;
[u1_prix,p1_prix,k1_prix,J1_prix,~,~] = DecompositionPrix(N1,A1,b1,C1,rho,eps,kmax,false);
[u1_alloc,omega1_alloc,k1_alloc,J1_alloc,~,~] = DecompositionQuantites(N1,A1,b1,C1,eps,kmax,false);
[u1_pred,omega1_pred,k1_pred,J1_pred,~,~] = DecompositionPrediction(N1,A1,b1,C1,eps,kmax,false);

%% Pour N = 3:
%   Solution: u* = (1;-1;0) et J(u*) = -1
N2 = 3;
[A2,b2,C2] = CreateInstance(N2);
rho = 0.1;
eps = 10^(-6);
kmax = 10000;
[u2_prix,p2_prix,k2_prix,J2_prix,~,~] = DecompositionPrix(N2,A2,b2,C2,rho,eps,kmax,false);
[u2_alloc,omega2_alloc,k2_alloc,J2_alloc,~,~] = DecompositionQuantites(N2,A2,b2,C2,eps,kmax,false);
[u2_pred,omega2_pred,k2_pred,J2_pred,~,~] = DecompositionPrediction(N2,A2,b2,C2,eps,kmax,false);

%% Pour N = 4:
%   Solution: u* = (1;-6/5;3/5;-1) et J(u*) = -19/10
N3 = 4;
[A3,b3,C3] = CreateInstance(N3);
rho = 0.1;
eps = 10^(-6);
kmax = 10000;
[u3_prix,p3_prix,k3_prix,J3_prix,~,~] = DecompositionPrix(N3,A3,b3,C3,rho,eps,kmax,false);
[u3_alloc,omega3_alloc,k3_alloc,J3_alloc,~,~] = DecompositionQuantites(N3,A3,b3,C3,eps,kmax,false);
[u3_pred,omega3_pred,k3_pred,J3_pred,~,~] = DecompositionPrediction(N3,A3,b3,C3,eps,kmax,false);

%% Pour N = 5:
%   Solution: u* = (1;-1.2;0.6;-1;0) et J(u*) = -1.9
N4 = 5;
[A4,b4,C4] = CreateInstance(N4);
rho = 0.1;
eps = 10^(-6);
kmax = 10000;

[u4_prix,p4_prix,k4_prix,J4_prix,~,~] = DecompositionPrix(N4,A4,b4,C4,rho,eps,kmax,false);
[u4_alloc,omega4_alloc,k4_alloc,J4_alloc,~,~] = DecompositionQuantites(N4,A4,b4,C4,eps,kmax,false);
[u4_pred,omega4_pred,k4_pred,J4_pred,~,~] = DecompositionPrediction(N4,A4,b4,C4,eps,kmax,false);






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 3
% On reprend les valeurs de A4, b4, C4 de l'encadré précédent
[~,p4_prix,~,~,~,U_prix] = DecompositionPrix(N4,A4,b4,C4,rho,eps,kmax,true);
[~,omega4_alloc,~,~,~,U_alloc] = DecompositionQuantites(N4,A4,b4,C4,eps,kmax,true);
[~,omega4_pred,~,~,~,U_pred] = DecompositionPrediction(N4,A4,b4,C4,eps,kmax,true);

%Pour tester à l'itération j :
j = 10;
disp(['Itération ',num2str(j),' :'])
disp(['Décomposition par les prix :',mat2str(Test_KKT(A4,b4,C4,0,p4_prix,0,0,0,U_prix(:,min(j,length(U_prix))),0.01))]);
disp(['Décomposition par les quantités :',mat2str(Test_KKT(A4,b4,C4,0,omega4_alloc,0,0,0,U_alloc(:,min(j,length(U_alloc))),0.01))]);
disp(['Décomposition par prédiction :',mat2str(Test_KKT(A4,b4,C4,0,omega4_pred,0,0,0,U_pred(:,min(j,length(U_pred))),0.01))]);


%On teste aux itérations j = 5,10,20,50,100,500 et 1000 des 3 algorithmes

for j = [5,10,20,50,100,500,1000]
    disp(['Itération ',num2str(j),' :'])
    disp(['Décomposition par les prix :',mat2str(Test_KKT(A4,b4,C4,0,p4_prix,0,0,0,U_prix(:,min(j,length(U_prix))),0.01))]);
    disp(['Décomposition par les quantités :',mat2str(Test_KKT(A4,b4,C4,0,omega4_alloc,0,0,0,U_alloc(:,min(j,length(U_alloc))),0.01))]);
    disp(['Décomposition par prédiction :',mat2str(Test_KKT(A4,b4,C4,0,omega4_pred,0,0,0,U_pred(:,min(j,length(U_pred))),0.01))]);
end







%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 4
rho = 0.1;
eps = 10^(-6);
kmax = 10000;

%%
plotNpetit = 0; %Si on veut visualiser la recherche de solution pour les 3 algos pour un petit exemple (N=5)
plotNbIt = 1; %Si on veut visualiser l'évolution du nombre d'itérations pour chaque algo pour N de 5 à 50
plotTmpsEx = 1; %Si on veut visualiser l'évolution du temps d'execution pour chaque algo pour N de 5 à 50
plotComplex = 1; %Si on veut visualiser l'erreur des algorithmes pour chaque algo pour N de 5 à 50

ComparMethode(rho,eps,kmax,plotNpetit,plotNbIt,plotTmpsEx,plotComplex)







%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 5
N5 = 200;
[A5,b5,C5] = CreateInstance(N5);
kmax = 10000;

%% Variation de epsilon
%%%%  Attention /!\ long à l'execution %%%%
rho = 0.1;
n=(6-2+1)*2;
eps_vec=zeros(n,1);
eps_vec(1:2:n)=10.^-(2:6);
eps_vec(2:2:n)=5*10.^-(3:7);
k5_prix = zeros(n,1); k5_alloc = zeros(n,1); k5_pred = zeros(n,1);
t_prix = zeros(n,1); t_alloc = zeros(n,1); t_pred = zeros(n,1);


for i=1:n
    [~,~,k5_prix(i),~,t_prix(i),~] = DecompositionPrix(N5,A5,b5,C5,rho,eps_vec(i),kmax,false);
    [~,~,k5_alloc(i),~,t_alloc(i),~] = DecompositionQuantites(N5,A5,b5,C5,rho,eps_vec(i),kmax,false);
    [~,~,k5_pred(i),~,t_pred(i),~] = DecompositionPrediction(N5,A5,b5,C5,eps_vec(i),kmax,false);
end

figure(1)
semilogx(eps_vec,[k5_prix,k5_alloc,k5_pred])
legend('Prix','Quantités','Prédiction')
xlabel('Précision (\epsilon)')
ylabel('Nombre d itérations (k)')

figure(2)
semilogx(eps_vec,[t_prix,t_alloc,t_pred])
legend('Prix','Quantités','Prédiction')
xlabel('Précision (\epsilon)')
ylabel('Temps d executon (en s)')

%% Variation de rho
%%%%  Attention /!\ long à l'execution %%%%
eps=10^-4;

rho_vec = 0.05:0.05:0.35;
n=length(rho_vec);
k5_prix = zeros(n,1); k5_alloc = zeros(n,1); %k5_pred = zeros(n,1);
t_prix = zeros(n,1); t_alloc = zeros(n,1); %t_pred = zeros(n,1);


for i=1:n
    [~,~,k5_prix(i),~,t_prix(i),~] = DecompositionPrix(N5,A5,b5,C5,rho_vec(i),eps,kmax,false);
    [~,~,k5_alloc(i),~,t_alloc(i),~] = DecompositionQuantites(N5,A5,b5,C5,rho_vec(i),eps,kmax,false);
    %[~,~,k5_pred(i),~,t_pred(i),~] = DecompositionPrediction(N5,A5,b5,C5,eps_vec(i),kmax,false);
end

figure(1)
plot(rho_vec,[k5_prix,k5_alloc])
legend('Prix','Quantités')
xlabel('Pas (\rho)')
ylabel('Nombre d itérations (k)')

figure(2)
plot(rho_vec,[t_prix,t_alloc])
legend('Prix','Quantités')
xlabel('Pas (\rho)')
ylabel('Temps d executon (en s)')

%%
%reste à tester de faire varier les paramètres alpha/bêta pour prédiction





