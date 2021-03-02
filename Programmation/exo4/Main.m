%SCENARIO 1 (Decomposition par les prix)

%Donnees du probleme:
Donnees_SmartGrids_Scenario1;

%Initialisation:
eps = 10^(-10);
kmax = 2000;
rho = 0.1;

%Algorithme de decomposition par les prix:
[P,~,k,J_opt] = DecompositionPrix(N,P0,a,b,Pmax,rho,eps,kmax);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 1');
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt)]);
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P(1)),' MW']);
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

%Solution exacte:
P_exact = [8.68; 2; 2; 2; 2; 2; 4.08; 2; 2; -9.98; -7.39; -8.98; -0.1; -0.11; -0.2];
J = @(P) sum(a.*(P - P0).^2 + b);

J_theo = J(P_exact);

%%
%SCENARIO 2 (Decomposition par les prix)

%Donnees du probleme:
Donnees_SmartGrids_Scenario2;

%Initialisation:
eps = 10^(-10);
kmax = 2000;
rho = 0.01;

%Algorithme de decomposition par les prix:
[P,p,k,J_opt] = DecompositionPrix(N,P0,a,b,Pmax,rho,eps,kmax);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 2');
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt)]);
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P(1)),' MW']);
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

%Solution exacte:
P_exact = [8.68; 2; 2; 2; 2; 2; 4.08; 2; 2; -9.98; -7.39; -8.98; -0.1; -0.11; -0.2];
J = @(P) sum(a.*(P - P0).^2 + b);

J_theo = J(P_exact);

%% 
%SCENARIO 1 (Decomposition par les quantités)

%Donnees du probleme:
Donnees_SmartGrids_Scenario1;

%Initialisation:
eps = 10^(-3);
kmax = 3000;

%Algorithme de decomposition par les prix:
[P,omega,k,J_opt] = DecompositionQuantites(N,P0,a,b,Pmax,eps,kmax);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 1');
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt)]);
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P(1)),' MW']);
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

%Solution exacte:
P_exact = [8.68; 2; 2; 2; 2; 2; 4.08; 2; 2; -9.98; -7.39; -8.98; -0.1; -0.11; -0.2];
J = @(P) sum(a.*(P - P0).^2 + b);

J_theo = J(P_exact);

%% 
%SCENARIO 2 (Decomposition par les quantités)

%Donnees du probleme:
Donnees_SmartGrids_Scenario2;

%Initialisation:
eps = 10^(-3);
kmax = 3000;

%Algorithme de decomposition par les prix:
[P,omega,k,J_opt] = DecompositionQuantites(N,P0,a,b,Pmax,eps,kmax);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 1');
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt)]);
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P(1)),' MW']);
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

%Solution exacte:
P_exact = [8.68; 2; 2; 2; 2; 2; 4.08; 2; 2; -9.98; -7.39; -8.98; -0.1; -0.11; -0.2];
J = @(P) sum(a.*(P - P0).^2 + b);

J_theo = J(P_exact);




