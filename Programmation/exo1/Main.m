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
set(objet,'facealpha',0.5);
hold off

title({'Ensemble des solutions','admissibles de J(u)'});
xlabel('u1');
ylabel('u2');
legend('u1 + 2u2 = 0','u2 = 0');

%% Question 2

%Test des decomposition par les prix, par les quantites et par prediction sur les exemples jouets:

%Pour N = 2: 
%   Solution: u* = (1;-1) et J(u*) = -1
N1 = 2;
[A1,b1,C1] = CreateInstance(N1);
rho = 0.1;
eps = 10^(-6);
kmax = 10000;
[u1_prix,p1_prix,k1_prix,J1_prix] = DecompositionPrix(N1,A1,b1,C1,rho,eps,kmax);
[u1_alloc,omega1_alloc,k1_alloc,J1_alloc] = DecompositionQuantites(N1,A1,b1,C1,eps,kmax);
[u1_pred,omega1_pred,k1_pred,J1_pred] = DecompositionPrediction(N1,A1,b1,C1,eps,kmax);

%Pour N = 3:
%   Solution: u* = (1;-1;0) et J(u*) = -1
N2 = 3;
[A2,b2,C2] = CreateInstance(N2);
rho = 0.1;
eps = 10^(-5);
kmax = 10000;
[u2_prix,p2_prix,k2_prix,J2_prix] = DecompositionPrix(N2,A2,b2,C2,rho,eps,kmax);
[u2_alloc,omega2_alloc,k2_alloc,J2_alloc] = DecompositionQuantites(N2,A2,b2,C2,eps,kmax);
[u2_pred,omega2_pred,k2_pred,J2_pred] = DecompositionPrediction(N2,A2,b2,C2,eps,kmax);

%Pour N = 4:
%   Solution: u* = (1;-6/5;3/5;-1) et J(u*) = -19/10
N3 = 4;
[A3,b3,C3] = CreateInstance(N3);
rho = 0.1;
eps = 10^(-6);
kmax = 10000;
[u3_prix,p3_prix,k3_prix,J3_prix] = DecompositionPrix(N3,A3,b3,C3,rho,eps,kmax);
[u3_alloc,omega3_alloc,k3_alloc,J3_alloc] = DecompositionQuantites(N3,A3,b3,C3,eps,kmax);
[u3_pred,omega3_pred,k3_pred,J3_pred] = DecompositionPrediction(N3,A3,b3,C3,eps,kmax);

%Pour N = 5:
%   Solution: u* = (1;-1.2;0.6;-1;0) et J(u*) = -1.9
N4 = 5;
[A4,b4,C4] = CreateInstance(N4);
rho = 0.1;
eps = 10^(-6);
kmax = 10000;
[u4_prix,p4_prix,k4_prix,J4_prix] = DecompositionPrix(N4,A4,b4,C4,rho,eps,kmax);
[u4_alloc,omega4_alloc,k4_alloc,J4_alloc] = DecompositionQuantites(N4,A4,b4,C4,eps,kmax);
[u4_pred,omega4_pred,k4_pred,J4_pred] = DecompositionPrediction(N4,A4,b4,C4,eps,kmax);

%% Question 3

Test_prix1 = Test_KKT(A1,b1,C1,0,p1_prix,0,0,0,u1_prix);
Test_prix2 = Test_KKT(A2,b2,C2,0,p2_prix,0,0,0,u2_prix);
Test_prix3 = Test_KKT(A3,b3,C3,0,p3_prix,0,0,0,u3_prix);
Test_prix4 = Test_KKT(A4,b4,C4,0,p4_prix,0,0,0,u4_prix);

