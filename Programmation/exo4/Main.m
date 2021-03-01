%Pour ajouter les algorithmes 
addpath('..\Algorithme');

%Donnees du probleme:
Donnees_SmartGrids;

%Initialisation:
eps = 10^(-8);
kmax = 2000;
p_ini = 0;
rho = 0.05;

%Algorithme de decomposition par les prix:
[P,p,k,J] = DecompositionPrix(N,P0,a,b,Pmax,rho,p_ini,eps,kmax);

%Affichage solution finale:
disp(['Valeur optimale: ',num2str(J)]);
disp('Puissances des agents:');
disp(['- Centrale_�_charbon: ',num2str(P(1)),' MW']);
disp(['- Eolienne 1: ',num2str(P(2)),' MW']);
disp(['- Eolienne 2: ',num2str(P(3)),' MW']);
disp(['- Eolienne 3: ',num2str(P(4)),' MW']);
disp(['- Eolienne 4: ',num2str(P(5)),' MW']);
disp(['- Eolienne 5: ',num2str(P(6)),' MW']);
disp(['- Barrage: ',num2str(P(7)),' MW']);
disp(['- Panneau photovoltaique 1: ',num2str(P(8)),' MW']);
disp(['- Panneau photovoltaique 2: ',num2str(P(9)),' MW']);
disp(['- DataCenter: ',num2str(P(10)),' MW']);
disp(['- Logement (x7500): ',num2str(P(11)),' MW']);
disp(['- Usine: ',num2str(P(12)),' MW']);
disp(['- Tramway 1: ',num2str(P(13)),' MW']);
disp(['- Tramway 2: ',num2str(P(14)),' MW']);
disp(['- Hopital: ',num2str(P(15)),' MW']);

