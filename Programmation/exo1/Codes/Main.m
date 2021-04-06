%EXERCICE 1:

%Antoine GICQUEL
%Cesar ARNOULD
%Hugo TESSIER
%Victor BERTRET

%% Question 1
clc;

%DONNEES GENERALES:

%Limites du domaine:
x_min = -10;
x_max = 10;
y_min = -10;
y_max = 10;

%Raffinement du domaine:
pas = 0.1;

%Fonction objectif à minimiser:
J = @(u1,u2) 0.5*(u1.^2 + u2.^2) - u1 + u2;

%Equations des droites delimitant le domaine:
f1 = @(x) - 1/2*x; f2 = @(x) x*0;
x1 = x_min:pas:x_max; y1 = f1(x1);
x2 = x_min:pas:x_max; y2 = f2(x2); 

%---------

%REPRESENTATION EN 2D:

%Positionnement de la figure:
figure(1);
set(gcf,'Position',[40 90 600 400]);

%Coloration de l'ensemble admissible:
ind = @(u1,u2) 1*((u1 + 2*u2 <= 0) & (u2 <= 0));
x_value = x_min:pas:x_max; y_value = y_min:pas:y_max;
[X,Y] = meshgrid(x_value,y_value);
z = ind(X,Y);
objet1_1 = fill([x_min x_min 0 x_max x_max],[y_min 0 0 f1(x_max) y_min],[0.73 0.82 0.88],'FaceAlpha',0.5);
hold on;

%Affichage des droites delimitant le domaine:
objet2_1 = plot(x1,y1,'r');
hold on;
objet3_1 = plot(x2,y2,'b');
hold on;

%Affichage de la solution exacte:
[u_solveur,~] = ResolutionExact(2);
objet4_1 = plot(u_solveur(1),u_solveur(2),'Color',[0 0 .5],'Marker','*');
hold on;

%Affichage des lignes de niveau de la fonction:
contour(x_value,y_value,J(X,Y),15,':');
colormap(winter(1000));
hold off

%Elements graphiques:
title({'Ensemble des solutions','admissibles de J(u)','(avec les lignes de niveau de J)'});
xlabel('u1');
ylabel('u2');
Lgd1 = legend([objet1_1,objet2_1,objet3_1,objet4_1],{'Ensemble admissible','u1 + 2u2 = 0','u2 = 0','Solution optimale'},'Location','EastOutside');
Lgd1.FontSize = 8;

%-----------------------

%REPRESENTATION EN 3D:

%Positionnement de la figure:
figure(2);
set(gcf,'Position',[650 90 600 400]);

%Coloration de l'ensemble admissible:
ind = @(u1,u2) 1*((u1 + 2*u2 <= 0) & (u2 <= 0));
x_value = x_min:pas:x_max; y_value = y_min:pas:y_max;
[X,Y] = meshgrid(x_value,y_value);
objet1_2 = fill([x_min x_min 0 x_max x_max],[y_min 0 0 f1(x_max) y_min],[0.73 0.82 0.88],'FaceAlpha',0.5);
hold on;

%Affichage de la fonction objectif en 3D:
Z = J(X,Y); Z_complementaire = J(X,Y);
mask = logical(ind(X,Y));
Z(~mask) = NaN; Z_complementaire(mask) = NaN;
objet2_2 = surf(X,Y,Z,'EdgeColor','none','FaceAlpha',0.8);
objet3_2 = surf(X,Y,Z_complementaire,'EdgeColor','none','FaceAlpha',0.2);
colormap(winter(1000));
caxis([-60 60]);
view(60,50);
hold on;

%Affichage des droites delimitant le domaine:
objet4_2 = plot(x1,y1,'r');
hold on;
objet5_2 = plot(x2,y2,'b');
hold on;

%Affichage de la solution exacte:
[u_solveur,~] = ResolutionExact(2);
objet6_2 = plot(u_solveur(1),u_solveur(2),'Color',[0 0 .5],'Marker','*');
hold on;

%Affichage des lignes de niveau de la fonction:
Z = J(X,Y);
contour(x_value,y_value,J(X,Y),15,':');
hold off

%Elements graphiques:
title({'Ensemble des solutions','admissibles de J(u)','(avec la surface 3D de J)'});
xlabel('u1');
ylabel('u2');
Lgd2 = legend([objet1_2,objet4_2,objet5_2,objet6_2,objet3_2],{'Ensemble admissible','u1 + 2u2 = 0','u2 = 0','Solution optimale','Surface (Admissible) de J'},'Location','EastOutside');
Lgd2.FontSize = 8;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 2

%IMPLEMENTATION DES ALGORITHMES DE DECOMPOSTION PAR PRIX, PAR QUANTITES, PAR PREDICTION:
% N = 2: 
clc;

N = 2;
[A,b,C] = CreateInstance(N);

%PARAMETRES GENERAUX:
rho = 0.1; %Pas pour la coordination
eps = 10^(-6); %Tolerance globale
kmax = 4000; %Nombre d'iterations maximal
PrintIt = true; %Affichage de chaque iteration dans la console
bigU = true; %Recuperation de toutes les solutions a chaque iteration

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax,...
                    'PrintIt',PrintIt,...
                    'bigU',bigU);
                
%PARAMETRES SOUS-PROBLEMES:
rho_sp_uzawa = 0.01; %Pas des sous-problemes
eps_sp = 10^(-4); %Tolerance des sous-problemes
kmax_sp = 1000; %Nombre d'iterations maximal pour les sous-problèmes

parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                                  'eps_sp',eps_sp,...
                                  'kmax_sp',kmax_sp);

%DECOMPOSITION PAR PRIX:
t = tic;
[u_prix,~,k_prix,J_prix] = DecompositionPrix(N,A,b,C,parametres,parametres_sousproblemes);
t_prix = toc(t);
%DECOMPOSITION PAR QUANTITES:
t = tic;
[u_alloc,~,k_alloc,J_alloc] = DecompositionQuantites(N,A,b,C,parametres,parametres_sousproblemes);
t_alloc = toc(t);
%DECOMPOSITION PAR PREDICTION:
t = tic;
[u_pred,~,k_pred,J_pred] = DecompositionPrediction(N,A,b,C,parametres,parametres_sousproblemes);
t_pred = toc(t);
[u_solveur,J_solveur] = ResolutionExact(N);

%AFFICHAGE:
disp(' ');
FonctionAffichage = {@(x) disp(['Solution = ( ', mat2str(round(x',3)),' )'])...,
                     @(J) disp(['Objectif = ', mat2str(J) ])...,
                     @(it) disp(['Nombres d''iterations = ', mat2str(it)])...,
                     @(t) disp(['Temps = ', mat2str(t) ,' s'])};
                 
disp('RESULTATS SOLVEUR:');
FonctionAffichage{1}(u_solveur(:,end));
FonctionAffichage{2}(J_solveur);
disp(' ');
disp('RESULTATS DECOMPOSITION PRIX:');
FonctionAffichage{1}(u_prix(:,end));
FonctionAffichage{2}(J_prix);
FonctionAffichage{3}(k_prix);
FonctionAffichage{4}(t_prix);
disp(' ');
disp('RESULTATS DECOMPOSITION QUANTITES:');
FonctionAffichage{1}(u_alloc(:,end));
FonctionAffichage{2}(J_alloc);
FonctionAffichage{3}(k_alloc);
FonctionAffichage{4}(t_alloc);
disp(' ');
disp('RESULTATS DECOMPOSITION PREDICTION:');
FonctionAffichage{1}(u_pred(:,end));
FonctionAffichage{2}(J_pred);
FonctionAffichage{3}(k_pred);
FonctionAffichage{4}(t_pred);

%% 

%IMPLEMENTATION DES ALGORITHMES DE DECOMPOSTION PAR PRIX, PAR QUANTITES, PAR PREDICTION:
% N = 3: 
clc;

N = 3;
[A,b,C] = CreateInstance(N);

%PARAMETRES GENERAUX:
rho = 0.1; %Pas pour la coordination
eps = 10^(-6); %Tolérance globale
kmax = 4000; %Nombre d'iterations maximal
PrintIt = true; %Affichage de chaque iteration dans la console
bigU = true; %Recuperation de toutes les solutions a chaque iteration

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax,...
                    'PrintIt',PrintIt,...
                    'bigU',bigU);
                
%PARAMETRES SOUS-PROBLEMES:
rho_sp_uzawa = 0.01; %Pas des sous-problemes
eps_sp = 10^(-4); %Tolerance des sous-problemes
kmax_sp = 1000; %Nombre d'iterations maximal pour les sous-problèmes

parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                                  'eps_sp',eps_sp,...
                                  'kmax_sp',kmax_sp);

%DECOMPOSITION PAR PRIX:
t = tic;
[u_prix,~,k_prix,J_prix] = DecompositionPrix(N,A,b,C,parametres,parametres_sousproblemes);
t_prix = toc(t);
%DECOMPOSITION PAR QUANTITES:
t = tic;
[u_alloc,~,k_alloc,J_alloc] = DecompositionQuantites(N,A,b,C,parametres,parametres_sousproblemes);
t_alloc = toc(t);
%DECOMPOSITION PAR PREDICTION:
t = tic;
[u_pred,~,k_pred,J_pred] = DecompositionPrediction(N,A,b,C,parametres,parametres_sousproblemes);
t_pred = toc(t);
[u_solveur,J_solveur] = ResolutionExact(N);

%AFFICHAGE:
disp(' ');
FonctionAffichage = {@(x) disp(['Solution = ( ', mat2str(round(x',3)),' )'])...,
                     @(J) disp(['Objectif = ', mat2str(J) ])...,
                     @(it) disp(['Nombres d''iterations = ', mat2str(it)])...,
                     @(t) disp(['Temps = ', mat2str(t) ,' s'])};
                 
disp('RESULTATS SOLVEUR:');
FonctionAffichage{1}(u_solveur(:,end));
FonctionAffichage{2}(J_solveur);
disp(' ');
disp('RESULTATS DECOMPOSITION PRIX:');
FonctionAffichage{1}(u_prix(:,end));
FonctionAffichage{2}(J_prix);
FonctionAffichage{3}(k_prix);
FonctionAffichage{4}(t_prix);
disp(' ');
disp('RESULTATS DECOMPOSITION QUANTITES:');
FonctionAffichage{1}(u_alloc(:,end));
FonctionAffichage{2}(J_alloc);
FonctionAffichage{3}(k_alloc);
FonctionAffichage{4}(t_alloc);
disp(' ');
disp('RESULTATS DECOMPOSITION PREDICTION:');
FonctionAffichage{1}(u_pred(:,end));
FonctionAffichage{2}(J_pred);
FonctionAffichage{3}(k_pred);
FonctionAffichage{4}(t_pred);

%%
%IMPLEMENTATION DES ALGORITHMES DE DECOMPOSTION PAR PRIX, PAR QUANTITES, PAR PREDICTION:
% N = 4: 
clc;

N = 4;
[A,b,C] = CreateInstance(N);

%PARAMETRES GENERAUX:
rho = 0.1; %Pas pour la coordination
eps = 10^(-6); %Tolérance globale
kmax = 4000; %Nombre d'iterations maximal
PrintIt = true; %Affichage de chaque iteration dans la console
bigU = true; %Recuperation de toutes les solutions a chaque iteration

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax,...
                    'PrintIt',PrintIt,...
                    'bigU',bigU);
                
%PARAMETRES SOUS-PROBLEMES:
rho_sp_uzawa = 0.01; %Pas des sous-problemes
eps_sp = 10^(-4); %Tolerance des sous-problemes
kmax_sp = 1000; %Nombre d'iterations maximal pour les sous-problèmes

parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                                  'eps_sp',eps_sp,...
                                  'kmax_sp',kmax_sp);

%DECOMPOSITION PAR PRIX:
t = tic;
[u_prix,~,k_prix,J_prix] = DecompositionPrix(N,A,b,C,parametres,parametres_sousproblemes);
t_prix = toc(t);
%DECOMPOSITION PAR QUANTITES:
t = tic;
[u_alloc,~,k_alloc,J_alloc] = DecompositionQuantites(N,A,b,C,parametres,parametres_sousproblemes);
t_alloc = toc(t);
%DECOMPOSITION PAR PREDICTION:
t = tic;
[u_pred,~,k_pred,J_pred] = DecompositionPrediction(N,A,b,C,parametres,parametres_sousproblemes);
t_pred = toc(t);
[u_solveur,J_solveur] = ResolutionExact(N);

%AFFICHAGE:
disp(' ');
FonctionAffichage = {@(x) disp(['Solution = ( ', mat2str(round(x',3)),' )'])...,
                     @(J) disp(['Objectif = ', mat2str(J) ])...,
                     @(it) disp(['Nombres d''iterations = ', mat2str(it)])...,
                     @(t) disp(['Temps = ', mat2str(t) ,' s'])};
                 
disp('RESULTATS SOLVEUR:');
FonctionAffichage{1}(u_solveur(:,end));
FonctionAffichage{2}(J_solveur);
disp(' ');
disp('RESULTATS DECOMPOSITION PRIX:');
FonctionAffichage{1}(u_prix(:,end));
FonctionAffichage{2}(J_prix);
FonctionAffichage{3}(k_prix);
FonctionAffichage{4}(t_prix);
disp(' ');
disp('RESULTATS DECOMPOSITION QUANTITES:');
FonctionAffichage{1}(u_alloc(:,end));
FonctionAffichage{2}(J_alloc);
FonctionAffichage{3}(k_alloc);
FonctionAffichage{4}(t_alloc);
disp(' ');
disp('RESULTATS DECOMPOSITION PREDICTION:');
FonctionAffichage{1}(u_pred(:,end));
FonctionAffichage{2}(J_pred);
FonctionAffichage{3}(k_pred);
FonctionAffichage{4}(t_pred);
%%
%IMPLEMENTATION DES ALGORITHMES DE DECOMPOSTION PAR PRIX, PAR QUANTITES, PAR PREDICTION:
% N = 5: 
clc;

N = 5;
[A,b,C] = CreateInstance(N);

%PARAMETRES GENERAUX:
rho = 0.1; %Pas pour la coordination
eps = 10^(-6); %Tolérance globale
kmax = 4000; %Nombre d'iterations maximal
PrintIt = true; %Affichage de chaque iteration dans la console
bigU = true; %Recuperation de toutes les solutions a chaque iteration

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax,...
                    'PrintIt',PrintIt,...
                    'bigU',bigU);
                
%PARAMETRES SOUS-PROBLEMES:
rho_sp_uzawa = 0.01; %Pas des sous-problemes
eps_sp = 10^(-4); %Tolerance des sous-problemes
kmax_sp = 1000; %Nombre d'iterations maximal pour les sous-problèmes

parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                                  'eps_sp',eps_sp,...
                                  'kmax_sp',kmax_sp);

%DECOMPOSITION PAR PRIX:
t = tic;
[u_prix,~,k_prix,J_prix] = DecompositionPrix(N,A,b,C,parametres,parametres_sousproblemes);
t_prix = toc(t);
%DECOMPOSITION PAR QUANTITES:
t = tic;
[u_alloc,~,k_alloc,J_alloc] = DecompositionQuantites(N,A,b,C,parametres,parametres_sousproblemes);
t_alloc = toc(t);
%DECOMPOSITION PAR PREDICTION:
t = tic;
[u_pred,~,k_pred,J_pred] = DecompositionPrediction(N,A,b,C,parametres,parametres_sousproblemes);
t_pred = toc(t);
[u_solveur,J_solveur] = ResolutionExact(N);

%AFFICHAGE:
disp(' ');
FonctionAffichage = {@(x) disp(['Solution = ( ', mat2str(round(x',3)),' )'])...,
                     @(J) disp(['Objectif = ', mat2str(J) ])...,
                     @(it) disp(['Nombres d''iterations = ', mat2str(it)])...,
                     @(t) disp(['Temps = ', mat2str(t) ,' s'])};
                 
disp('RESULTATS SOLVEUR:');
FonctionAffichage{1}(u_solveur(:,end));
FonctionAffichage{2}(J_solveur);
disp(' ');
disp('RESULTATS DECOMPOSITION PRIX:');
FonctionAffichage{1}(u_prix(:,end));
FonctionAffichage{2}(J_prix);
FonctionAffichage{3}(k_prix);
FonctionAffichage{4}(t_prix);
disp(' ');
disp('RESULTATS DECOMPOSITION QUANTITES:');
FonctionAffichage{1}(u_alloc(:,end));
FonctionAffichage{2}(J_alloc);
FonctionAffichage{3}(k_alloc);
FonctionAffichage{4}(t_alloc);
disp(' ');
disp('RESULTATS DECOMPOSITION PREDICTION:');
FonctionAffichage{1}(u_pred(:,end));
FonctionAffichage{2}(J_pred);
FonctionAffichage{3}(k_pred);
FonctionAffichage{4}(t_pred);

%%
%IMPLEMENTATION DES ALGORITHMES DE DECOMPOSTION PAR PRIX, PAR QUANTITES, PAR PREDICTION:
% N = 10: 
clc;

N = 10;
[A,b,C] = CreateInstance(N);

%PARAMETRES GENERAUX:
rho = 0.1; %Pas pour la coordination
eps = 10^(-6); %Tolérance globale
kmax = 4000; %Nombre d'iterations maximal
PrintIt = true; %Affichage de chaque iteration dans la console
bigU = true; %Recuperation de toutes les solutions a chaque iteration

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax,...
                    'PrintIt',PrintIt,...
                    'bigU',bigU);
                
%PARAMETRES SOUS-PROBLEMES:
rho_sp_uzawa = 0.01; %Pas des sous-problemes
eps_sp = 10^(-4); %Tolerance des sous-problemes
kmax_sp = 1000; %Nombre d'iterations maximal pour les sous-problèmes

parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                                  'eps_sp',eps_sp,...
                                  'kmax_sp',kmax_sp);

%DECOMPOSITION PAR PRIX:
t = tic;
[u_prix,~,k_prix,J_prix] = DecompositionPrix(N,A,b,C,parametres,parametres_sousproblemes);
t_prix = toc(t);
%DECOMPOSITION PAR QUANTITES:
t = tic;
[u_alloc,~,k_alloc,J_alloc] = DecompositionQuantites(N,A,b,C,parametres,parametres_sousproblemes);
t_alloc = toc(t);
%DECOMPOSITION PAR PREDICTION:
t = tic;
[u_pred,~,k_pred,J_pred] = DecompositionPrediction(N,A,b,C,parametres,parametres_sousproblemes);
t_pred = toc(t);
[u_solveur,J_solveur] = ResolutionExact(N);

%AFFICHAGE:
disp(' ');
FonctionAffichage = {@(x) disp(['Solution = ( ', mat2str(round(x',3)),' )'])...,
                     @(J) disp(['Objectif = ', mat2str(J) ])...,
                     @(it) disp(['Nombres d''iterations = ', mat2str(it)])...,
                     @(t) disp(['Temps = ', mat2str(t) ,' s'])};
                 
disp('RESULTATS SOLVEUR:');
FonctionAffichage{1}(u_solveur(:,end));
FonctionAffichage{2}(J_solveur);
disp(' ');
disp('RESULTATS DECOMPOSITION PRIX:');
FonctionAffichage{1}(u_prix(:,end));
FonctionAffichage{2}(J_prix);
FonctionAffichage{3}(k_prix);
FonctionAffichage{4}(t_prix);
disp(' ');
disp('RESULTATS DECOMPOSITION QUANTITES:');
FonctionAffichage{1}(u_alloc(:,end));
FonctionAffichage{2}(J_alloc);
FonctionAffichage{3}(k_alloc);
FonctionAffichage{4}(t_alloc);
disp(' ');
disp('RESULTATS DECOMPOSITION PREDICTION:');
FonctionAffichage{1}(u_pred(:,end));
FonctionAffichage{2}(J_pred);
FonctionAffichage{3}(k_pred);
FonctionAffichage{4}(t_pred);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 3
clc;

%Pour une valeur fixée de N, on vérifie d'abord que KKT est bien verifié à
%l'issue de la convergence de l'algorithme

N =10;
[A,b,C] = CreateInstance(N);

%PARAMETRES GENERAUX:
rho = 0.1; %Pas pour la coordination
eps = 10^(-6); %Tolérance globale
kmax = 4000; %Nombre d'iterations maximal
PrintIt = false; %Affichage de chaque iteration dans la console
bigU = false; %Recuperation de toutes les solutions a chaque iteration

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax,...
                    'PrintIt',PrintIt,...
                    'bigU',bigU);
                
%PARAMETRES SOUS-PROBLEMES:
rho_sp_uzawa = 0.01; %Pas des sous-problemes
eps_sp = 10^(-4); %Tolerance des sous-problemes
kmax_sp = 1000; %Nombre d'iterations maximale pour les sous-problemes

parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                                  'eps_sp',eps_sp,...
                                  'kmax_sp',kmax_sp);
                              
CheckKKT = {@(x) disp(['- Contraintes d''inégalité ? : ',mat2str(x)])...,
            @(x) disp(['- Multiplicateurs Mu positifs ? : ',mat2str(x)])...,
            @(x) disp(['- Multiplicateurs Mu * Contraintes d''inégalité nul ? : ',mat2str(x)])...,
            @(x) disp(['- Contraintes égalité ? : ',mat2str(x)])...,
            @(x) disp(['- Gradient Lagrangien nul ? : ',mat2str(x)])};

%TEST DES CONDITIONS DE KKT A LA FIN DE CHAQUE ALGORITHME:
disp('TEST KKT À L''ISSUE DE CHAQUE ALGORITHME :');

%Decomposition par les prix:
[u_prix,Mu_prix] = DecompositionPrix(N,A,b,C,parametres,parametres_sousproblemes);
[OK_prix,cond1_prix,cond2_prix,cond3_prix,cond4_prix,cond5_prix] = Test_KKT(A,b,C,zeros(N,1),Mu_prix,0,0,0,u_prix,0.001);
disp(['Conditions KKT verifiées pour la décomposition par prix ? : ',mat2str(OK_prix)]);
CheckKKT{1}(cond1_prix);
CheckKKT{2}(cond2_prix);
CheckKKT{3}(cond3_prix);
CheckKKT{4}(cond4_prix);
CheckKKT{5}(cond5_prix);
disp(' ');

%Decomposition par quantites:
[u_alloc,Mu_alloc] = DecompositionQuantites(N,A,b,C,parametres,parametres_sousproblemes);
[OK_alloc,cond1_alloc,cond2_alloc,cond3_alloc,cond4_alloc,cond5_alloc] = Test_KKT(A,b,C,zeros(N,1),Mu_alloc,0,0,0,u_alloc,0.001);
disp(['Conditions KKT verifiées pour la décomposition par quantités ? : ',mat2str(OK_alloc)]);
CheckKKT{1}(cond1_alloc);
CheckKKT{2}(cond2_alloc);
CheckKKT{3}(cond3_alloc);
CheckKKT{4}(cond4_alloc);
CheckKKT{5}(cond5_alloc);
disp(' ');

%Decomposition par prediction:
[u_pred,Mu_pred,~,~,~,~] = DecompositionPrediction(N,A,b,C,parametres,parametres_sousproblemes);
[OK_pred,cond1_pred,cond2_pred,cond3_pred,cond4_pred,cond5_pred] = Test_KKT(A,b,C,zeros(N,1),Mu_pred,0,0,0,u_pred,0.001);
disp(['Conditions KKT verifiées pour la décomposition par prédiction ? : ',mat2str(OK_pred)]);
CheckKKT{1}(cond1_pred);
CheckKKT{2}(cond2_pred);
CheckKKT{3}(cond3_pred);
CheckKKT{4}(cond4_pred);
CheckKKT{5}(cond5_pred);
disp(' ');
%%
clc;
N = 8;
[A,b,C] = CreateInstance(N);

%PARAMETRES GENERAUX:
rho = 0.1; %Pas pour la coordination
eps = 10^(-6); %Tolérance globale
kmax = 4000; %Nombre d'iterations maximal
PrintIt = false; %Affichage de chaque iteration dans la console
bigU = true; %Recuperation de toutes les solutions a chaque iteration
bigMu = true; %Recuperation des multiplicateurs a chaque iteration

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax,...
                    'PrintIt',PrintIt,...
                    'bigU',bigU,...
                    'bigMu',bigMu);
                
%PARAMETRES SOUS-PROBLEMES:
rho_sp_uzawa = 0.01; %Pas des sous-problemes
eps_sp = 10^(-4); %Tolerance des sous-problemes
kmax_sp = 1000; %Nombre d'iterations maximal pour les sous-problèmes

parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                                  'eps_sp',eps_sp,...
                                  'kmax_sp',kmax_sp);

%TEST DES CONDITIONS DE KKT A UNE ITERATION DONNEE:

%Decomposition par les prix:
[u_prix,Mu_prix,~,~,~,U_prix_final,Mu_prix_final] = DecompositionPrix(N,A,b,C,parametres,parametres_sousproblemes);

%Decomposition par quantites:
[u_alloc,Mu_alloc,~,~,~,U_alloc_final,Mu_alloc_final] = DecompositionQuantites(N,A,b,C,parametres,parametres_sousproblemes);

%Decomposition par prediction:
[u_pred,Mu_pred,~,~,~,U_pred_final,Mu_pred_final] = DecompositionPrediction(N,A,b,C,parametres,parametres_sousproblemes);


CheckKKT = {@(x) disp(['- Contraintes d''inégalité ? : ',mat2str(x)])...,
            @(x) disp(['- Multiplicateurs Mu positifs ? : ',mat2str(x)])...,
            @(x) disp(['- Multiplicateurs Mu * Contraintes d''inégalité nul ? : ',mat2str(x)])...,
            @(x) disp(['- Contraintes égalité ? : ',mat2str(x)])...,
            @(x) disp(['- Gradient Lagrangien nul ? : ',mat2str(x)])};
                 
j = 10;
disp(['Itération ',mat2str(j),' :']);
%Decomposition par les prix:
[OK_prixj,cond1_prixj,cond2_prixj,cond3_prixj,cond4_prixj,cond5_prixj] = Test_KKT(A,b,C,zeros(N,1),Mu_prix_final(:,j),0,0,0,U_prix_final(:,j),0.001);
disp(['Conditions KKT verifiées pour la decomposition par prix ? : ',mat2str(OK_prixj)]);
CheckKKT{1}(cond1_prixj);
CheckKKT{2}(cond2_prixj);
CheckKKT{3}(cond3_prixj);
CheckKKT{4}(cond4_prixj);
CheckKKT{5}(cond5_prixj);
disp(' ');
%Decomposition par les quantités:
[OK_allocj,cond1_allocj,cond2_allocj,cond3_allocj,cond4_allocj,cond5_allocj] = Test_KKT(A,b,C,zeros(N,1),Mu_alloc_final(:,j),0,0,0,U_alloc_final(:,j),0.001);
disp(['Conditions KKT verifiées pour la decomposition par allocation ? : ',mat2str(OK_allocj)]);
CheckKKT{1}(cond1_allocj);
CheckKKT{2}(cond2_allocj);
CheckKKT{3}(cond3_allocj);
CheckKKT{4}(cond4_allocj);
CheckKKT{5}(cond5_allocj);
disp(' ');
%Decomposition par prediction:
[OK_predj,cond1_predj,cond2_predj,cond3_predj,cond4_predj,cond5_predj] = Test_KKT(A,b,C,zeros(N,1),Mu_pred_final(:,j),0,0,0,U_pred_final(:,j),0.001);
disp(['Conditions KKT verifiées pour la decomposition par prédiction ? : ',mat2str(OK_predj)]);
CheckKKT{1}(cond1_predj);
CheckKKT{2}(cond2_predj);
CheckKKT{3}(cond3_predj);
CheckKKT{4}(cond4_predj);
CheckKKT{5}(cond5_predj);
disp(' ');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 4
clc;

%COMPARAISON DES ALGORITHMES (TEMPS D'EXECUTION, NOMBRE D'ITERATIONS)

%OPTIONS:
plotNpetit = 1; %Si on veut visualiser la recherche de solution pour les 3 algos pour un petit exemple (N=5)
plotNbIt = 1; %Si on veut visualiser l'évolution du nombre d'iterations pour chaque algo pour N de 10 a 50
plotTmpsEx = 1; %Si on veut visualiser l'evolution du temps d'execution pour chaque algo pour N de 10 a 50
plotComplex = 1; %Si on veut visualiser l'erreur des algorithmes pour chaque algo pour N de 10 a 50

%PARAMETRES GENERAUX:
rho = 0.1; %Pas de la coordination
eps = 10^(-6); %Tolerance generale
kmax = 10000; %Nombre d'iterations maximal
PrintIt = false; %Affichage de chaque iteration dans la console
bigU = true; %Recuperation de toutes les solutions a chaque iteration

parametres = struct('rho', rho, ...
                    'eps', eps, ...
                    'kmax', kmax,...
                    'PrintIt',PrintIt,...
                    'bigU',bigU);
                
%PARAMETRES SOUS-PROBLEMES:
rho_sp_uzawa = 0.01; %Pas des sous-problemes
eps_sp = 10^(-4); %Tolerance des sous-problemes
kmax_sp = 1000; %Nombre d'iterations maximal pour les sous-problemes

parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                                  'eps_sp',eps_sp,...
                                  'kmax_sp',kmax_sp);

%COMPARAISON DES METHODES:
ComparMethode(parametres,parametres_sousproblemes,plotNpetit,plotNbIt,plotTmpsEx,plotComplex)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 5

%% Variation de epsilon
%%%%  Attention /!\ long a l'execution %%%%
clc;

%CREATION DE L'INSTANCE:
N = 200;
%N = 251;
[A,b,C] = CreateInstance(N);

%PARAMETRES GENERAUX:
rho = 0.1; %Pas de la coordination
kmax = 10000; %Nombre d'iterations maximal
PrintIt = false; %Affichage de chaque iteration dans la console
bigU = true; %Recuperation de toutes les solutions a chaque iteration

%PARAMETRES SOUS-PROBLEMES:
rho_sp_uzawa = 0.01; %Pas des sous-problemes
eps_sp = 10^(-4); %Tolerance des sous-problemes
kmax_sp = 1000; %Nombre d'iterations maximal pour les sous-problèmes

%VALEURS DES EPSILONS:
n=(6-2+1)*2;
eps_vec=zeros(n,1);
eps_vec(1:2:n)=10.^-(2:6);
eps_vec(2:2:n)=5*10.^-(3:7);

%VALEURS A MESURER:
k5_prix = zeros(n,1); k5_alloc = zeros(n,1); k5_pred = zeros(n,1);
t_prix = zeros(n,1); t_alloc = zeros(n,1); t_pred = zeros(n,1);

disp('Exécution en cours...');

for i=1:n
    disp(['n = ',mat2str(i),' | epsilon = ',mat2str(eps_vec(i))]);
    
    %Initialisation generale :
    eps = eps_vec(i);
    parametres = struct('rho', rho, ...
                        'eps', eps, ...
                        'kmax', kmax,...
                        'PrintIt',PrintIt,...
                        'bigU',bigU);

    %Initialisation des sous-problemes :
    parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                                      'eps_sp',eps_sp,...
                                      'kmax_sp',kmax_sp);
    
    %Decomposition par les prix :
    [~,~,k5_prix(i),~,t_prix(i),~] = DecompositionPrix(N,A,b,C,parametres,parametres_sousproblemes);
    disp(' - Décomposition Prix OK');
    %Decomposition par quantites :
    [~,~,k5_alloc(i),~,t_alloc(i),~] = DecompositionQuantites(N,A,b,C,parametres,parametres_sousproblemes);
    disp(' - Décomposition Quantités OK');
    %Decomposition par prediction :
    [~,~,k5_pred(i),~,t_pred(i),~] = DecompositionPrediction(N,A,b,C,parametres,parametres_sousproblemes);
    disp(' - Décomposition Prédiction OK');
end

%AFFICHAGE EVOLUTION DU NOMBRE D'ITERATIONS
figure(1)
semilogx(eps_vec,[k5_prix,k5_alloc,k5_pred]);
title({'Evolution du nombre d''itérations','en fonction des valeurs de precision'});
legend('Décomposition par prix','Décomposition par quantités','Décomposition par prédiction');
xlabel('Précision (\epsilon)');
ylabel('Nombre d itérations (k)');

%AFFICHAGE EVOLUTION DU TEMPS D'EXECUTION:
figure(2)
semilogx(eps_vec,[t_prix,t_alloc,t_pred])
title({'Evolution du temps d''exécution','en fonction des valeurs de précision'});
legend('Décomposition par prix','Décomposition par quantités','Décomposition par prédiction');
xlabel('Précision (\epsilon)');
ylabel('Temps d''éxécuton (en s)');

%% Variation de rho
%%%%  Attention /!\ long a l'execution %%%%
clc;

%CREATION DE L'INSTANCE:
N = 200;
%N = 251;
[A,b,C] = CreateInstance(N);

%PARAMETRES GENERAUX:
eps = 10^(-4); %Tolerance generale
kmax = 10000; %Nombre d'iterations maximal
PrintIt = false; %Affichage de chaque iteration dans la console
bigU = true; %Recuperation de toutes les solutions a chaque iteration

%PARAMETRES SOUS-PROBLEMES:
rho_sp_uzawa = 0.01; %Pas des sous-problemes
eps_sp = 10^(-4); %Tolerance des sous-problemes
kmax_sp = 1000; %Nombre d'iterations maximal pour les sous-problemes

%VALEURS DES RHOS:
%rho_vec = 0.05:0.05:0.35;
rho_vec = [0.005 0.01 0.025 0.05 0.1 0.15 0.2 0.25 0.3 0.35];

%VALEURS A MESURER:
n=length(rho_vec);
k5_prix = zeros(n,1); k5_alloc = zeros(n,1); %k5_pred = zeros(n,1);
t_prix = zeros(n,1); t_alloc = zeros(n,1); %t_pred = zeros(n,1);

disp('Exécution en cours...');

for i=1:n
    disp(['n = ',mat2str(i),' | rho = ',mat2str(rho_vec(i))]);
    
    %Initialisation général:
    rho = rho_vec(i);
    parametres = struct('rho', rho, ...
                        'eps', eps, ...
                        'kmax', kmax,...
                        'PrintIt',PrintIt,...
                        'bigU',bigU);

    %Initialisation des sous-problemes:
    parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
                                      'eps_sp',eps_sp,...
                                      'kmax_sp',kmax_sp);
    
    %Decomposition par les prix:
    [~,~,k5_prix(i),~,t_prix(i),~] = DecompositionPrix(N,A,b,C,parametres,parametres_sousproblemes);
    disp(' - Décomposition Prix OK');
    %Decomposition par quantites:
    [~,~,k5_alloc(i),~,t_alloc(i),~] = DecompositionQuantites(N,A,b,C,parametres,parametres_sousproblemes);
    disp(' - Décomposition Quantités OK');
end

%AFFICHAGE EVOLUTION DU NOMBRE D'ITERATIONS
figure(1)
plot(rho_vec,[k5_prix,k5_alloc]);
title({'Evolution du nombre d''itérations','en fonction des valeurs de pas'});
legend('Décomposition par prix','Décomposition par quantités');
xlabel('Pas (\rho)');
ylabel('Nombre d itérations (k)');

%AFFICHAGE EVOLUTION DU TEMPS D'EXECUTION:
figure(2)
plot(rho_vec,[t_prix,t_alloc])
title({'Evolution du temps d''exécution','en fonction des valeurs de pas'});
legend('Décomposition par prix','Décomposition par quantités');
xlabel('Pas (\rho)');
ylabel('Temps d''exécuton (en s)');

%% Variation des parametres de relaxation
%%%%  Attention /!\ long a l'execution %%%%
clc;

%CREATION DE L'INSTANCE:
N = 200;
%N = 251;
[A,b,C] = CreateInstance(N);

%PARAMETRES GENERAUX:
eps = 10^(-4); %Tolerance generale
kmax = 10000; %Nombre d'iterations maximal
PrintIt = false; %Affichage de chaque iteration dans la console
bigU = false; %Recuperation de toutes les solutions a chaque iteration

%PARAMETRES SOUS-PROBLEMES:
rho_sp_uzawa = 0.01; %Pas des sous-problemes
eps_sp = 10^(-4); %Tolerance des sous-problemes
kmax_sp = 1000; %Nombre d'iterations maximal pour les sous-problemes

%VALEURS DES PARAMETRES DE RELAXATION:
beta_vec = 0.1:0.1:0.9;
gamma_vec = 0.1:0.1:0.9;

%VALEURS A MESURER:
n=length(beta_vec);
m=length(gamma_vec);

k5_pred = zeros(n,m);
t5_pred = zeros(n,m);

disp('Exécution en cours...');

for i=1:n
    for j = 1:m
        disp(['i = ',mat2str(i),' | j = ',mat2str(j),' | bêta = ',mat2str(beta_vec(i)),' | gamma = ',mat2str(gamma_vec(j))]);
        
        %Initialisation general:
        beta = beta_vec(i);
        gamma = gamma_vec(j);
        parametres = struct('eps', eps, ...
            'kmax', kmax,...
            'PrintIt',PrintIt,...
            'bigU',bigU,...
            'beta',beta,...
            'gamma',gamma);
        
        %Initialisation des sous-problemes:
        parametres_sousproblemes = struct('rho_sp_uzawa',rho_sp_uzawa,...
            'eps_sp',eps_sp,...
            'kmax_sp',kmax_sp);
        
        %Decomposition par les prix:
        [~,~,k5_pred(i,j),~,t5_pred(i,j),~] = DecompositionPrediction(N,A,b,C,parametres,parametres_sousproblemes);
        disp(' - Décomposition Prédiction OK');
    end
end

% AFFICHAGE EVOLUTION DU NOMBRE D'ITERATIONS
figure(1)
pcolor(beta_vec,gamma_vec,k5_pred);
colorbar;
title({'Evolution du nombre d''itérations','en fonction des valeurs de gamma et de bêta','(Décomposition par prédiction)'});
xlabel('gamma');
ylabel('bêta');

%AFFICHAGE EVOLUTION DU TEMPS D'EXECUTION:
figure(2)
pcolor(beta_vec,gamma_vec,t5_pred);
colorbar;
title({'Evolution du temps d''exécution','en fonction des valeurs de gamma et de bêta','(Décomposition par prédiction)'});
xlabel('gamma');
ylabel('bêta');
