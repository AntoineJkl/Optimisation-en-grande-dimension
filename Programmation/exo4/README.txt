/Codes : répertoire des codes correspondant à l'exercice 4

	- Main.m: Fichier principal. Il se compose de 4 cellules indépendantes. Dans chacune de ces cellules, on exécute l'algorithme correspondant, on affiche la solution dans la console ainsi que l'évolution de la solution et des valeurs du multiplicateur en fonction des itérations de l'algorithme:
		- La première cellule correspond à l'exécution de l'algorithme de décomposition par les prix pour le scénario 1
		- La deuxième cellule correspond à l'exécution de l'algorithme de décomposition par les prix pour le scénario 2
		- La troisième cellule correspond à l'exécution de l'algorithme de décomposition par les quantités pour le scénario 1
		- La quatrième cellule correspond à l'exécution de l'algorithme de décomposition par les quantités pour le scénario 2
	
	- Donnees_SmartGrids_Scenario1.m: Fichier contenant toutes les données de l'instance du scénario 1 (a,b,P0,Pmax pour chaque agent).

	- Donnees_SmartGrids_Scenario2.m: Fichier contenant toutes les données de l'instance du scénario 2 (a,b,P0,Pmax pour chaque agent).

	- DecompositionPrix.m: Implémentation de l'algorithme de décomposition par les prix. La fonction prend en entrée les données de l'instance ainsi qu'un jeu de paramètres et renvoie en sortie la solution obtenue (à chaque iteration), le nombre d'itérations effectuées, le temps d'exécution, la valeur objective en la solution obtenue J(P*) et les valeurs du multiplicateur à chaque itération.

	- DecompositionQuantites.m: Implémentation de l'algorithme de décomposition par les quantités. La fonction prend en entrée les données de l'instance ainsi qu'un jeu de paramètres et renvoie en sortie la solution obtenue (à chaque iteration), le nombre d'itérations effectuées, le temps d'exécution, la valeur objective en la solution obtenue J(P*) et les valeurs du multiplicateur à chaque itération.

	- ResolutionExact.m: Fonction permettant de résoudre une instance données en paramètre à l'aide du solveur de MATLAB. Elle prend en entrée les données de l'instance à traiter et renvoie la solution, la valeur objective et la valeur du multiplicateur. Cette fonction nous sert de référence.

	- AffichageSolution.m: Fonction permettant d'afficher les valeurs des solutions en fonction des itérations lors de l'exécution de l'algorithme de décomposition par les prix ou par les quantités.

	- AffichageMultiplicateur.m: Fonction permettant d'afficher les valeurs du multiplicateur de Lagrange associé à la contrainte d'égalité du problème en fonction des itérations lors de l'exécution de l'algorithme de décomposition par les prix ou par les quantités.



/Figures : répertoire des figures présentées dans le rapport. 