%NOMBRE DE SOUS-PROBLEMES:
N = 15;

%DONNEES SUR LES AGENTS:
P0 = zeros(N,1);
a = zeros(N,1);
b = zeros(N,1);
Pmax = zeros(N,1);
noms = ['Centrale_à_charbon',...
        'Eolienne_1',... 
        'Eolienne_2',...
        'Eolienne_3',...
        'Eolienne_4',...
        'Eolienne_5',...
        'Barrage',...
        'Panneau_photovoltaique_1',...
        'Panneau_photovoltaique_2',...
        'DataCenter',...
        'Logement_(x7500)',...
        'Usine',...
        'Tramway_1',...
        'Tramway_2',...
        'Hopital'];

%Centrale a charbon:
P0(1) = 8.5;
a(1) = 10;
b(1) = 1000;
Pmax(1) = 10;

%Eolienne 1:
P0(2) = 1.9;
a(2) = 1;
b(2) = 10;
Pmax(2) = 2;

%Eolienne 2:
P0(3) = 1.9;
a(3) = 1;
b(3) = 10;
Pmax(3) = 2;

%Eolienne 3:
P0(4) = 1.9;
a(4) = 1;
b(4) = 10;
Pmax(4) = 2;

%Eolienne 4:
P0(5) = 1.9;
a(5) = 1;
b(5) = 10;
Pmax(5) = 2;

%Eolienne 5:
P0(6) = 1.9;
a(6) = 1;
b(6) = 10;
Pmax(6) = 2;

%Barrage:
P0(7) = 3;
a(7) = 0.1;
b(7) = 100;
Pmax(7) = 16;

%Panneau photovoltaique 1:
P0(8) = 1.9;
a(8) = 0.9;
b(8) = 11;
Pmax(8) = 2;

%Panneau photovoltaique 2:
P0(9) = 1.9;
a(9) = 0.9;
b(9) = 11;
Pmax(9) = 2;

%DataCenter:
P0(10) = -10;
a(10) = 5;
b(10) = 20;
Pmax(10) = 0;

%Logement (7500):
P0(11) = -7.5;
a(11) = 1;
b(11) = 10;
Pmax(11) = 0;

%Usine:
P0(12) = -9;
a(12) = 5;
b(12) = 8;
Pmax(12) = 0;

%Tramway 1:
P0(13) = -0.12;
a(13) = 6;
b(13) = 9;
Pmax(13) = 0;

%Tramway 2:
P0(14) = -0.12;
a(14) = 6;
b(14) = 9;
Pmax(14) = 0;

%Hopital:
P0(15) = -0.2;
a(15) = 200;
b(15) = 30;
Pmax(15) = 0;

