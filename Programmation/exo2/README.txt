- Main.m: Fichier principal. Il se compose de 20 cellules indépendantes. Voici le résumé de chaque cellule : 
	- La 1ère cellule correspond à l'exécution des algorithmes de décomposition sur l'exemple Jouet avec Arrow ou Uzawa
	- La 2ème cellule correspond à l'exécution des algorithmes de décomposition par quantité sur l'exemple Jouet avec les formules explicites
	- La 3ème cellule correspond à l'exécution des algorithmes de décomposition sur l'exemple Jouet avec la méthode des points intérieurs
	- Les 4,5 et 6èmes cellules correspondent aux affichages des résultats des cellules précédentes
	- Les 7 et 8èmes cellules correspondent à l'exécution des algorithmes de décomposition par prédiction sur l'exemple Jouet sans relaxation
	- Les 9 et 10èmes cellules correspondent à l'exécution de l'algorithme de décomposition par prédiction seq sur l'exemple Jouet avec relaxation sur l'allocation
	- Les 11 et 12èmes cellules correspondent à l'étude de l'influence de beta et gamma sur les algorithmes de décomposition par prédiction
	- Les 13 et 14èmes cellules correspondent à l'étude du nombre d'itérations et du temps d'exécution des algorithmes de décomposition par prédiction en fonction de N
	- Les 15 et 16èmes cellules correspondent à l'étude du nombre d'iterations et du temps d'execution en fonction du pas de la descente
	- Les 17 et 18èmes cellules correspondent à l'étude du nombre d'iterations et du temps d'execution en fonction de N
	- Les 19 et 20èmes cellules correspondent à l'étude du nombre d'iterations et du temps d'execution en fonction de la précision

-CreationInstance.m : Implémentation d'une fonction créant des instances aléatoires en fonction des valeurs maximales et minimales souhaitées

-Juv.m : Fonction objective associée aux sous problèmes. Elle prend en paramètre la matrice A, le vecteur b et var la valeur avec laquelle on veut calculer la valeur de la fonction objective.


/decomposition : répertoire des décompositions

	- DecompositionPrix.m: Implémentation de l'algorithme de décomposition par les prix. La fonction prend en entrée un jeu de paramètres et renvoie en sortie la solution obtenue , le prix optimal , le nombre d'itérations effectuées, la valeur objective finale et l'ensemble des solutions u à chaque itération.

	- DecompositionQuantites.m: Implémentation de l'algorithme de décomposition par les quantités. La fonction prend en entrée un jeu de paramètres et renvoie en sortie la solution obtenue , les allocations optimaux , le nombre d'itérations effectuées, la valeur objective finale et l'ensemble des solutions u à chaque itération.

	- DecompositionQuantites2.m: Implémentation de l'algorithme de décomposition par les quantités 2. La fonction prend en entrée un jeu de paramètres et renvoie en sortie la solution obtenue , les allocations optimaux , le nombre d'itérations effectuées, la valeur objective finale et l'ensemble des solutions u à chaque itération.

	- DecompositionPredictionPar.m: Implémentation de l'algorithme de décomposition par prédiction parallèle. La fonction prend en entrée un jeu de paramètres et renvoie en sortie la solution obtenue , les allocations optimaux , le prix optimal , le nombre d'itérations effectuées, la valeur objective finale et l'ensemble des solutions u à chaque itération.	
	
	- DecompositionPredictionSeq.m: Implémentation de l'algorithme de décomposition par prédiction séquentielle. La fonction prend en entrée un jeu de paramètres et renvoie en sortie la solution obtenue , les allocations optimaux , le prix optimal , le nombre d'itérations effectuées, la valeur objective finale et l'ensemble des solutions u à chaque itération.

	- DecompositionPredictionSeq2.m: Implémentation de l'algorithme de décomposition par prédiction séquentielle 2. La fonction prend en entrée un jeu de paramètres et renvoie en sortie la solution obtenue , les allocations optimaux , le prix optimal , le nombre d'itérations effectuées, la valeur objective finale et l'ensemble des solutions u à chaque itération.


/Figures : répertoire des figures présentées dans le rapport. 