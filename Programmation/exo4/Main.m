% SCENARIO 1 (Decomposition par les prix)

%Donnees du probleme:
Donnees_SmartGrids_Scenario1;

%Initialisation:
eps = 10^(-4);
kmax = 10000;
rho = 0.1;

%Algorithme de decomposition par les prix:
[P,~,k,J_opt] = DecompositionPrix(N,P0,a,b,Pmax,rho,eps,kmax);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 1 - Décomposition par les prix');
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt)]);
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P(1,end)),' MW']);
disp(['- Eolienne 1: ',num2str(P(2,end)),' MW']);
disp(['- Eolienne 2: ',num2str(P(3,end)),' MW']);
disp(['- Eolienne 3: ',num2str(P(4,end)),' MW']);
disp(['- Eolienne 4: ',num2str(P(5,end)),' MW']);
disp(['- Eolienne 5: ',num2str(P(6,end)),' MW']);
disp(['- Barrage: ',num2str(P(7,end)),' MW']);
disp(['- Panneau photovoltaique 1: ',num2str(P(8,end)),' MW']);
disp(['- Panneau photovoltaique 2: ',num2str(P(9,end)),' MW']);
disp(['- DataCenter: ',num2str(P(10,end)),' MW']);
disp(['- Logement (x7500): ',num2str(P(11,end)),' MW']);
disp(['- Usine: ',num2str(P(12,end)),' MW']);
disp(['- Tramway 1: ',num2str(P(13,end)),' MW']);
disp(['- Tramway 2: ',num2str(P(14,end)),' MW']);
disp(['- Hopital: ',num2str(P(15,end)),' MW']);

%Solution exacte:
P_exact = [8.68; 2; 2; 2; 2; 2; 4.08; 2; 2; -9.98; -7.39; -8.98; -0.1; -0.11; -0.2];
J = @(P) sum(a.*(P - P0).^2 + b);

J_theo = J(P_exact);

%%
%SCENARIO 2 (Decomposition par les prix)

%Donnees du probleme:
Donnees_SmartGrids_Scenario2;

%Initialisation:
eps = 10^(-4);
kmax = 10000;
rho = 0.01;

%Algorithme de decomposition par les prix:
[P,p,k,J_opt] = DecompositionPrix(N,P0,a,b,Pmax,rho,eps,kmax);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 2 - Décomposition par les prix');
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt)]);
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P(1,end)),' MW']);
disp(['- Eolienne 1: ',num2str(P(2,end)),' MW']);
disp(['- Eolienne 2: ',num2str(P(3,end)),' MW']);
disp(['- Eolienne 3: ',num2str(P(4,end)),' MW']);
disp(['- Eolienne 4: ',num2str(P(5,end)),' MW']);
disp(['- Eolienne 5: ',num2str(P(6,end)),' MW']);
disp(['- Barrage: ',num2str(P(7,end)),' MW']);
disp(['- Panneau photovoltaique 1: ',num2str(P(8,end)),' MW']);
disp(['- Panneau photovoltaique 2: ',num2str(P(9,end)),' MW']);
disp(['- DataCenter: ',num2str(P(10,end)),' MW']);
disp(['- Logement (x7500): ',num2str(P(11,end)),' MW']);
disp(['- Usine: ',num2str(P(12,end)),' MW']);
disp(['- Tramway 1: ',num2str(P(13,end)),' MW']);
disp(['- Tramway 2: ',num2str(P(14,end)),' MW']);
disp(['- Hopital: ',num2str(P(15,end)),' MW']);

%Solution exacte:
P_exact = [8.68; 2; 2; 2; 2; 2; 4.08; 2; 2; -9.98; -7.39; -8.98; -0.10; -0.11; -0.2];
J = @(P) sum(a.*(P - P0).^2 + b);

J_theo = J(P_exact);

%% 
%SCENARIO 1 (Decomposition par les quantités)

%Donnees du probleme:
Donnees_SmartGrids_Scenario1;

%Initialisation:
eps = 10^(-4);
kmax = 10000;

%Algorithme de decomposition par les prix:
[P,omega,k,J_opt] = DecompositionQuantites(N,P0,a,b,Pmax,eps,kmax);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 1 - Décomposition par les quantités');
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt)]);
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P(1,end)),' MW']);
disp(['- Eolienne 1: ',num2str(P(2,end)),' MW']);
disp(['- Eolienne 2: ',num2str(P(3,end)),' MW']);
disp(['- Eolienne 3: ',num2str(P(4,end)),' MW']);
disp(['- Eolienne 4: ',num2str(P(5,end)),' MW']);
disp(['- Eolienne 5: ',num2str(P(6,end)),' MW']);
disp(['- Barrage: ',num2str(P(7,end)),' MW']);
disp(['- Panneau photovoltaique 1: ',num2str(P(8,end)),' MW']);
disp(['- Panneau photovoltaique 2: ',num2str(P(9,end)),' MW']);
disp(['- DataCenter: ',num2str(P(10,end)),' MW']);
disp(['- Logement (x7500): ',num2str(P(11,end)),' MW']);
disp(['- Usine: ',num2str(P(12,end)),' MW']);
disp(['- Tramway 1: ',num2str(P(13,end)),' MW']);
disp(['- Tramway 2: ',num2str(P(14,end)),' MW']);
disp(['- Hopital: ',num2str(P(15,end)),' MW']);

%Solution exacte:
P_exact = [8.68; 2; 2; 2; 2; 2; 4.08; 2; 2; -9.98; -7.39; -8.98; -0.1; -0.11; -0.2];
J = @(P) sum(a.*(P - P0).^2 + b);

J_theo = J(P_exact);

%% 
%SCENARIO 2 (Decomposition par les quantités)

%Donnees du probleme:
Donnees_SmartGrids_Scenario2;

%Initialisation:
eps = 10^(-4);
kmax = 10000;

%Algorithme de decomposition par les prix:
[P,omega,k,J_opt] = DecompositionQuantites(N,P0,a,b,Pmax,eps,kmax);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 2 - Décomposition par les quantités');
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt)]);
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P(1,end)),' MW']);
disp(['- Eolienne 1: ',num2str(P(2,end)),' MW']);
disp(['- Eolienne 2: ',num2str(P(3,end)),' MW']);
disp(['- Eolienne 3: ',num2str(P(4,end)),' MW']);
disp(['- Eolienne 4: ',num2str(P(5,end)),' MW']);
disp(['- Eolienne 5: ',num2str(P(6,end)),' MW']);
disp(['- Barrage: ',num2str(P(7,end)),' MW']);
disp(['- Panneau photovoltaique 1: ',num2str(P(8,end)),' MW']);
disp(['- Panneau photovoltaique 2: ',num2str(P(9,end)),' MW']);
disp(['- DataCenter: ',num2str(P(10,end)),' MW']);
disp(['- Logement (x7500): ',num2str(P(11,end)),' MW']);
disp(['- Usine: ',num2str(P(12,end)),' MW']);
disp(['- Tramway 1: ',num2str(P(13,end)),' MW']);
disp(['- Tramway 2: ',num2str(P(14,end)),' MW']);
disp(['- Hopital: ',num2str(P(15,end)),' MW']);

%Solution exacte:
P_exact = [9.93; 0; 0; 0; 0; 0; 11.59; 2; 2; -9.82; -6.64; -8.82; 0; 0.03; -0.2];
J = @(P) sum(a.*(P - P0).^2 + b);

J_theo = J(P_exact);




