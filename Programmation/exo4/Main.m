% SCENARIO 1 (Decomposition par les prix)
clc;

%Donnees du probleme:
Donnees_SmartGrids_Scenario1;

%Initialisation du problème général:
rho = 0.1;
eps = 10^(-5);
kmax = 10000;

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax);

%Initialisation des sous-problèmes:
choix = 1; %1 = Uzawa, 0 = Arrow
rho_sp_uzawa = 0.1;
rho_sp_arrow1 = 0.001;
rho_sp_arrow2 = 0.001;
eps_sp = 10^(-4);
kmax_sp = 10000;

parametres_sousproblemes = struct('choix',choix,...
                    'rho_sp_uzawa',rho_sp_uzawa,...
                    'rho_sp_arrow1',rho_sp_arrow1,...
                    'rho_sp_arrow2',rho_sp_arrow2,...
                    'eps_sp',eps_sp,...
                    'kmax_sp',kmax_sp);
                    
%Algorithme de decomposition par les prix:
[P,~,k,J_opt,time] = DecompositionPrix(N,P0,a,b,Pmax,parametres,parametres_sousproblemes);
P_final = P(:,end);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 1 - Décomposition par les prix');
disp(['Temps execution (s): ',num2str(time)]);
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt)]);
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P_final(1)),' MW']);
disp(['- Eolienne 1: ',num2str(P_final(2)),' MW']);
disp(['- Eolienne 2: ',num2str(P_final(3)),' MW']);
disp(['- Eolienne 3: ',num2str(P_final(4)),' MW']);
disp(['- Eolienne 4: ',num2str(P_final(5)),' MW']);
disp(['- Eolienne 5: ',num2str(P_final(6)),' MW']);
disp(['- Barrage: ',num2str(P_final(7)),' MW']);
disp(['- Panneau photovoltaique 1: ',num2str(P_final(8)),' MW']);
disp(['- Panneau photovoltaique 2: ',num2str(P_final(9)),' MW']);
disp(['- DataCenter: ',num2str(P_final(10)),' MW']);
disp(['- Logement (x7500): ',num2str(P_final(11)),' MW']);
disp(['- Usine: ',num2str(P_final(12)),' MW']);
disp(['- Tramway 1: ',num2str(P_final(13)),' MW']);
disp(['- Tramway 2: ',num2str(P_final(14)),' MW']);
disp(['- Hopital: ',num2str(P_final(15)),' MW']);

%Solution exacte:
[P_exact,J_exact] = ResolutionExact(N,P0,a,b,Pmax);

%Affichage des solutions:
AffichageSolution(P,P_exact);

%%
%SCENARIO 2 (Decomposition par les prix)
clc;

%Donnees du probleme:
Donnees_SmartGrids_Scenario2;

%Initialisation du problème général:
rho = 0.1;
eps = 10^(-4);
kmax = 10000;

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax);

%Initialisation des sous-problèmes:
choix = 1; %1 = Uzawa, 0 = Arrow
rho_sp_uzawa = 0.1;
rho_sp_arrow1 = 0.001;
rho_sp_arrow2 = 0.001;
eps_sp = 10^(-4);
kmax_sp = 10000;

parametres_sousproblemes = struct('choix',choix,...
                    'rho_sp_uzawa',rho_sp_uzawa,...
                    'rho_sp_arrow1',rho_sp_arrow1,...
                    'rho_sp_arrow2',rho_sp_arrow2,...
                    'eps_sp',eps_sp,...
                    'kmax_sp',kmax_sp);

%Algorithme de decomposition par les prix:
[P,~,k,J_opt,time] = DecompositionPrix(N,P0,a,b,Pmax,parametres,parametres_sousproblemes);
P_final = P(:,end);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 2 - Décomposition par les prix');
disp(['Temps execution (s): ',num2str(time)]);
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt)]);
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P_final(1)),' MW']);
disp(['- Eolienne 1: ',num2str(P_final(2)),' MW']);
disp(['- Eolienne 2: ',num2str(P_final(3)),' MW']);
disp(['- Eolienne 3: ',num2str(P_final(4)),' MW']);
disp(['- Eolienne 4: ',num2str(P_final(5)),' MW']);
disp(['- Eolienne 5: ',num2str(P_final(6)),' MW']);
disp(['- Barrage: ',num2str(P_final(7)),' MW']);
disp(['- Panneau photovoltaique 1: ',num2str(P_final(8)),' MW']);
disp(['- Panneau photovoltaique 2: ',num2str(P_final(9)),' MW']);
disp(['- DataCenter: ',num2str(P_final(10)),' MW']);
disp(['- Logement (x7500): ',num2str(P_final(11)),' MW']);
disp(['- Usine: ',num2str(P_final(12)),' MW']);
disp(['- Tramway 1: ',num2str(P_final(13)),' MW']);
disp(['- Tramway 2: ',num2str(P_final(14)),' MW']);
disp(['- Hopital: ',num2str(P_final(15)),' MW']);

%Solution exacte:
[P_exact,J_exact] = ResolutionExact(N,P0,a,b,Pmax);

%Affichage des solutions:
AffichageSolution(P,P_exact);

%% 
%SCENARIO 1 (Decomposition par les quantités)
clc;

%Donnees du probleme:
Donnees_SmartGrids_Scenario1;

%Initialisation du problème général:
rho = 0.1;
eps = 10^(-4);
kmax = 10000;

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax);

%Initialisation des sous-problèmes:
choix = 1; %1 = Uzawa, 0 = Arrow
rho_sp_uzawa = 0.05;
rho_sp_arrow1 = 0.001;
rho_sp_arrow2 = 0.001;
eps_sp = 10^(-3);
kmax_sp = 5000;

parametres_sousproblemes = struct('choix',choix,...
                    'rho_sp_uzawa',rho_sp_uzawa,...
                    'rho_sp_arrow1',rho_sp_arrow1,...
                    'rho_sp_arrow2',rho_sp_arrow2,...
                    'eps_sp',eps_sp,...
                    'kmax_sp',kmax_sp);

%Algorithme de decomposition par les prix:
[P,~,k,J_opt,time] = DecompositionQuantites(N,P0,a,b,Pmax,parametres,parametres_sousproblemes);
P_final = P(:,end);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 1 - Décomposition par quantités');
disp(['Temps execution (s): ',num2str(time)]);
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt)]);
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P_final(1)),' MW']);
disp(['- Eolienne 1: ',num2str(P_final(2)),' MW']);
disp(['- Eolienne 2: ',num2str(P_final(3)),' MW']);
disp(['- Eolienne 3: ',num2str(P_final(4)),' MW']);
disp(['- Eolienne 4: ',num2str(P_final(5)),' MW']);
disp(['- Eolienne 5: ',num2str(P_final(6)),' MW']);
disp(['- Barrage: ',num2str(P_final(7)),' MW']);
disp(['- Panneau photovoltaique 1: ',num2str(P_final(8)),' MW']);
disp(['- Panneau photovoltaique 2: ',num2str(P_final(9)),' MW']);
disp(['- DataCenter: ',num2str(P_final(10)),' MW']);
disp(['- Logement (x7500): ',num2str(P_final(11)),' MW']);
disp(['- Usine: ',num2str(P_final(12)),' MW']);
disp(['- Tramway 1: ',num2str(P_final(13)),' MW']);
disp(['- Tramway 2: ',num2str(P_final(14)),' MW']);
disp(['- Hopital: ',num2str(P_final(15)),' MW']);

%Solution exacte:
[P_exact,J_exact] = ResolutionExact(N,P0,a,b,Pmax);

%Affichage des solutions:
AffichageSolution(P,P_exact);

%% 
%SCENARIO 2 (Decomposition par les quantités)
clc;

%Donnees du probleme:
Donnees_SmartGrids_Scenario2;

%Initialisation du problème général:
rho = 0.1;
eps = 10^(-4);
kmax = 10000;

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax);

%Initialisation des sous-problèmes:
choix = 1; %1 = Uzawa, 0 = Arrow
rho_sp_uzawa = 0.005;
rho_sp_arrow1 = 0.001;
rho_sp_arrow2 = 0.001;
eps_sp = 10^(-3);
kmax_sp = 5000;

parametres_sousproblemes = struct('choix',choix,...
                    'rho_sp_uzawa',rho_sp_uzawa,...
                    'rho_sp_arrow1',rho_sp_arrow1,...
                    'rho_sp_arrow2',rho_sp_arrow2,...
                    'eps_sp',eps_sp,...
                    'kmax_sp',kmax_sp);

%Algorithme de decomposition par les prix:
[P,~,k,J_opt,time] = DecompositionQuantites(N,P0,a,b,Pmax,parametres,parametres_sousproblemes);
P_final = P(:,end);

%Affichage solution finale:
disp(' ');
disp('SCENARIO 2 - Décomposition par quantités');
disp(['Temps execution (s): ',num2str(time)]);
disp(['Nombre d iterations: ',num2str(k)]);
disp(['Valeur optimale: ',num2str(J_opt)]);
disp('Puissances des agents:');
disp(['- Centrale à charbon: ',num2str(P_final(1)),' MW']);
disp(['- Eolienne 1: ',num2str(P_final(2)),' MW']);
disp(['- Eolienne 2: ',num2str(P_final(3)),' MW']);
disp(['- Eolienne 3: ',num2str(P_final(4)),' MW']);
disp(['- Eolienne 4: ',num2str(P_final(5)),' MW']);
disp(['- Eolienne 5: ',num2str(P_final(6)),' MW']);
disp(['- Barrage: ',num2str(P_final(7)),' MW']);
disp(['- Panneau photovoltaique 1: ',num2str(P_final(8)),' MW']);
disp(['- Panneau photovoltaique 2: ',num2str(P_final(9)),' MW']);
disp(['- DataCenter: ',num2str(P_final(10)),' MW']);
disp(['- Logement (x7500): ',num2str(P_final(11)),' MW']);
disp(['- Usine: ',num2str(P_final(12)),' MW']);
disp(['- Tramway 1: ',num2str(P_final(13)),' MW']);
disp(['- Tramway 2: ',num2str(P_final(14)),' MW']);
disp(['- Hopital: ',num2str(P_final(15)),' MW']);

%Solution exacte:
[P_exact,J_exact] = ResolutionExact(N,P0,a,b,Pmax);

%Affichage des solutions:
AffichageSolution(P,P_exact);



