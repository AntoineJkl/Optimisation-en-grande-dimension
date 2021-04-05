% SCENARIO 1 (Decomposition par les prix)
clc;

%Donnees du probleme:
Donnees_SmartGrids_Scenario1;

%Initialisation du problème général:
rho = 0.05;
eps = 10^(-4);
kmax = 10000;

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax);

%Initialisation des sous-problèmes:
rho_sp_uzawa = 0.1;
eps_sp = 10^(-4);
kmax_sp = 5000;

parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                    'eps_sp',eps_sp,...
                    'kmax_sp',kmax_sp);
                    
%Algorithme de decomposition par les prix:
[P,Lambdas,k,J_opt,time] = DecompositionPrix(N,P0,a,b,Pmax,parametres,parametres_sousproblemes);
P_final = P(:,end);

%Solution exacte:
[P_exact,J_exact,Lambda_exact] = ResolutionExact(N,P0,a,b,Pmax);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 1 - Décomposition par les prix');
disp(['Temps execution (s): ',num2str(time)]);
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt),'               ( Valeur optimale solveur: ',num2str(J_exact),' )']);
disp(' ');
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P_final(1)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(1)),' MW )']);
disp(['- Eolienne 1: ',num2str(P_final(2)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(2)),' MW )']);
disp(['- Eolienne 2: ',num2str(P_final(3)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(3)),' MW )']);
disp(['- Eolienne 3: ',num2str(P_final(4)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(4)),' MW )']);
disp(['- Eolienne 4: ',num2str(P_final(5)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(5)),' MW )']);
disp(['- Eolienne 5: ',num2str(P_final(6)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(6)),' MW )']);
disp(['- Barrage: ',num2str(P_final(7)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(7)),' MW )']);
disp(['- Panneau photovoltaique 1: ',num2str(P_final(8)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(8)),' MW )']);
disp(['- Panneau photovoltaique 2: ',num2str(P_final(9)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(9)),' MW )']);
disp(['- DataCenter: ',num2str(P_final(10)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(10)),' MW )']);
disp(['- Logement (x7500): ',num2str(P_final(11)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(11)),' MW )']);
disp(['- Usine: ',num2str(P_final(12)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(12)),' MW )']);
disp(['- Tramway 1: ',num2str(P_final(13)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(13)),' MW )']);
disp(['- Tramway 2: ',num2str(P_final(14)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(14)),' MW )']);
disp(['- Hopital: ',num2str(P_final(15)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(15)),' MW )']);
disp(' ');
disp('Multiplicateur de la contrainte d''égalité: ');
disp(['- Lambda: ',num2str(Lambdas(end)),'               ( Valeur optimale solveur: ',num2str(Lambda_exact),' )']);


%Affichage multiplicateur:
titre1 = {'Convergence du multiplicateur','de la contrainte dans le scénario 1','(Décomposition par les prix)'};
AffichageMultiplicateur(Lambdas,Lambda_exact,titre1);

%Affichage des solutions:
titre2 = {'Convergence des solutions','pour chaque agent dans le scénario 1','(Décomposition par les prix)'};
AffichageSolution(P,P_exact,titre2);

%%
%SCENARIO 2 (Decomposition par les prix)
clc;

%Donnees du probleme:
Donnees_SmartGrids_Scenario2;

%Initialisation du problème général:
rho = 0.05;
eps = 10^(-4);
kmax = 10000;

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax);

%Initialisation des sous-problèmes:
rho_sp_uzawa = 0.1;
eps_sp = 10^(-4);
kmax_sp = 5000;

parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                    'eps_sp',eps_sp,...
                    'kmax_sp',kmax_sp);

%Algorithme de decomposition par les prix:
[P,Lambdas,k,J_opt,time] = DecompositionPrix(N,P0,a,b,Pmax,parametres,parametres_sousproblemes);
P_final = P(:,end);

%Solution exacte:
[P_exact,J_exact,Lambda_exact] = ResolutionExact(N,P0,a,b,Pmax);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 2 - Décomposition par les prix');
disp(['Temps execution (s): ',num2str(time)]);
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt),'               ( Valeur optimale solveur: ',num2str(J_exact),' )']);
disp(' ');
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P_final(1)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(1)),' MW )']);
disp(['- Eolienne 1: ',num2str(P_final(2)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(2)),' MW )']);
disp(['- Eolienne 2: ',num2str(P_final(3)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(3)),' MW )']);
disp(['- Eolienne 3: ',num2str(P_final(4)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(4)),' MW )']);
disp(['- Eolienne 4: ',num2str(P_final(5)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(5)),' MW )']);
disp(['- Eolienne 5: ',num2str(P_final(6)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(6)),' MW )']);
disp(['- Barrage: ',num2str(P_final(7)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(7)),' MW )']);
disp(['- Panneau photovoltaique 1: ',num2str(P_final(8)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(8)),' MW )']);
disp(['- Panneau photovoltaique 2: ',num2str(P_final(9)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(9)),' MW )']);
disp(['- DataCenter: ',num2str(P_final(10)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(10)),' MW )']);
disp(['- Logement (x7500): ',num2str(P_final(11)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(11)),' MW )']);
disp(['- Usine: ',num2str(P_final(12)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(12)),' MW )']);
disp(['- Tramway 1: ',num2str(P_final(13)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(13)),' MW )']);
disp(['- Tramway 2: ',num2str(P_final(14)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(14)),' MW )']);
disp(['- Hopital: ',num2str(P_final(15)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(15)),' MW )']);
disp(' ');
disp('Multiplicateur de la contrainte d''égalité: ');
disp(['- Lambda: ',num2str(Lambdas(end)),'               ( Valeur optimale solveur: ',num2str(Lambda_exact),' )']);

%Affichage multiplicateur:
titre1 = {'Convergence du multiplicateur','de la contrainte dans le scénario 2','(Décomposition par les prix)'};
AffichageMultiplicateur(Lambdas,Lambda_exact,titre1);

%Affichage des solutions:
titre2 = {'Convergence des solutions','pour chaque agent dans le scénario 2','(Décomposition par les prix)'};
AffichageSolution(P,P_exact,titre2);

%% 
%SCENARIO 1 (Decomposition par les quantités)
clc;

%Donnees du probleme:
Donnees_SmartGrids_Scenario1;

%Initialisation du problème général:
rho = 0.1;
eps = 10^(-4);
kmax = 50000;

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax);

%Initialisation des sous-problèmes:
rho_sp_uzawa = 0.01;
eps_sp = 10^(-4);
kmax_sp = 5000;

parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                    'eps_sp',eps_sp,...
                    'kmax_sp',kmax_sp);

%Algorithme de decomposition par les prix:
[P,Lambdas,k,J_opt,time] = DecompositionQuantites(N,P0,a,b,Pmax,parametres,parametres_sousproblemes);
P_final = P(:,end);

%Solution exacte:
[P_exact,J_exact,Lambda_exact] = ResolutionExact(N,P0,a,b,Pmax);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 1 - Décomposition par quantités');
disp(['Temps execution (s): ',num2str(time)]);
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt),'               ( Valeur optimale solveur: ',num2str(J_exact),' )']);
disp(' ');
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P_final(1)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(1)),' MW )']);
disp(['- Eolienne 1: ',num2str(P_final(2)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(2)),' MW )']);
disp(['- Eolienne 2: ',num2str(P_final(3)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(3)),' MW )']);
disp(['- Eolienne 3: ',num2str(P_final(4)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(4)),' MW )']);
disp(['- Eolienne 4: ',num2str(P_final(5)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(5)),' MW )']);
disp(['- Eolienne 5: ',num2str(P_final(6)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(6)),' MW )']);
disp(['- Barrage: ',num2str(P_final(7)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(7)),' MW )']);
disp(['- Panneau photovoltaique 1: ',num2str(P_final(8)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(8)),' MW )']);
disp(['- Panneau photovoltaique 2: ',num2str(P_final(9)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(9)),' MW )']);
disp(['- DataCenter: ',num2str(P_final(10)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(10)),' MW )']);
disp(['- Logement (x7500): ',num2str(P_final(11)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(11)),' MW )']);
disp(['- Usine: ',num2str(P_final(12)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(12)),' MW )']);
disp(['- Tramway 1: ',num2str(P_final(13)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(13)),' MW )']);
disp(['- Tramway 2: ',num2str(P_final(14)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(14)),' MW )']);
disp(['- Hopital: ',num2str(P_final(15)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(15)),' MW )']);
disp(' ');
disp('Multiplicateur de la contrainte d''égalité: ');
disp(['- Lambda: ',num2str(Lambdas(end)),'               ( Valeur optimale solveur: ',num2str(Lambda_exact),' )']);

%Affichage multiplicateur:
titre1 = {'Convergence du multiplicateur','de la contrainte dans le scénario 1','(Décomposition par les quantités)'};
AffichageMultiplicateur(Lambdas,Lambda_exact,titre1);

%Affichage des solutions:
titre2 = {'Convergence des solutions','pour chaque agent dans le scénario 1','(Décomposition par les quantités)'};
AffichageSolution(P,P_exact,titre2);

%% 
%SCENARIO 2 (Decomposition par les quantités)
clc;

%Donnees du probleme:
Donnees_SmartGrids_Scenario2;

%Initialisation du problème général:
rho = 0.01;
eps = 10^(-4);
kmax = 20000;

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax);

%Initialisation des sous-problèmes:
rho_sp_uzawa = 0.001;
eps_sp = 10^(-4);
kmax_sp = 5000;

parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                    'eps_sp',eps_sp,...
                    'kmax_sp',kmax_sp);

%Algorithme de decomposition par les prix:
[P,Lambdas,k,J_opt,time] = DecompositionQuantites(N,P0,a,b,Pmax,parametres,parametres_sousproblemes);
P_final = P(:,end);

%Solution exacte:
[P_exact,J_exact,Lambda_exact] = ResolutionExact(N,P0,a,b,Pmax);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 2 - Décomposition par quantités');
disp(['Temps execution (s): ',num2str(time)]);
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt),'               ( Valeur optimale solveur: ',num2str(J_exact),' )']);
disp(' ');
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P_final(1)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(1)),' MW )']);
disp(['- Eolienne 1: ',num2str(P_final(2)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(2)),' MW )']);
disp(['- Eolienne 2: ',num2str(P_final(3)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(3)),' MW )']);
disp(['- Eolienne 3: ',num2str(P_final(4)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(4)),' MW )']);
disp(['- Eolienne 4: ',num2str(P_final(5)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(5)),' MW )']);
disp(['- Eolienne 5: ',num2str(P_final(6)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(6)),' MW )']);
disp(['- Barrage: ',num2str(P_final(7)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(7)),' MW )']);
disp(['- Panneau photovoltaique 1: ',num2str(P_final(8)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(8)),' MW )']);
disp(['- Panneau photovoltaique 2: ',num2str(P_final(9)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(9)),' MW )']);
disp(['- DataCenter: ',num2str(P_final(10)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(10)),' MW )']);
disp(['- Logement (x7500): ',num2str(P_final(11)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(11)),' MW )']);
disp(['- Usine: ',num2str(P_final(12)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(12)),' MW )']);
disp(['- Tramway 1: ',num2str(P_final(13)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(13)),' MW )']);
disp(['- Tramway 2: ',num2str(P_final(14)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(14)),' MW )']);
disp(['- Hopital: ',num2str(P_final(15)),' MW','               ( Valeur optimale solveur: ',num2str(P_exact(15)),' MW )']);
disp(' ');
disp('Multiplicateur de la contrainte d''égalité: ');
disp(['- Lambda: ',num2str(Lambdas(end)),'               ( Valeur optimale solveur: ',num2str(Lambda_exact),' )']);

%Affichage multiplicateur:
titre1 = {'Convergence du multiplicateur','de la contrainte dans le scénario 2','(Décomposition par les quantités)'};
AffichageMultiplicateur(Lambdas,Lambda_exact,titre1);

%Affichage des solutions:
titre2 = {'Convergence des solutions','pour chaque agent dans le scénario 2','(Décomposition par les quantités)'};
AffichageSolution(P,P_exact,titre2);
