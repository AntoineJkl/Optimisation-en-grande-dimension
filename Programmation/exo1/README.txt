/Codes : répertoire des codes correspondant à l'exercice 1

	- Main.m: Fichier principal. Il se compose de plusieurs cellules indépendantes correspondant à chaque questions du sujet. Les questions sont indiquées comme titre de cellule. Si rien n'est indiqué, il s'agit de la même question que la précédente mais avec une variation de N ou des paramètres. EN executant chaque cellule, et celon les cas, on affiche dans la console l'évolution de l'exécution de chaque algorithme et/ou un graphe est tracé pour résumer les données et leur évolution :
		- La première cellule correspond à la première question : on représente en 2D et 3D l'ensemble admissible du problème dans le cas N=2
		- Les cellules 2 à 5 correspondent à la seconde question : on teste nos 3 algorithmes dans les cas N = 3, 4, 5 et 10 respectivement
		- Les cellules 6 et 7 correspondent à la troisème question : on teste si les conditions KKT sont vérifiées après convergence et à l'itération j=10 respectivement
		- La huitième cellule correspond à la quatrième question : on compare les 3 algorithmes en traçant au choix un petit exemple (N=5) d'atteinte des solutions, les variations de temps d'exécution et du nombre d'itérations, l'erreur commise par les algorithmes par rapoort à la solution exacte
		- Les cellules 9 à 11 correspondent à la cinquième question : on  fait varier dans le cas N=200 (ou N=251 en dé-commentant la ligne correspondante) les paramètres epsilon, rho et les paramètres de relaxation respectivement

	- ComparMethode.m: Fonction permettant la comparaison des méthodes décrite pour la question 4. Elle prend en entrée deux jeux de paramètres ainsi que des booléens pour indiquer si on calcule et affiche les différents items au choix.

	- CreateInstance.m: Fichier contenant toutes les données de l'instance pour notre problème. Elle crée l'instance correspondant à la taille N du problème.

	- DecompositionPrix.m: Implémentation de l'algorithme de décomposition par les prix. La fonction prend en entrée les données de l'instance ainsi qu'un jeu de paramètres et renvoie en sortie la solution obtenue (à chaque iteration), le nombre d'itérations effectuées, le temps d'exécution, la valeur objective en la solution obtenue J(P*), l'ensemble des solutions calculées à chaque itération et les valeurs du multiplicateur à chaque itération.

	- DecompositionQuantites.m: Implémentation de l'algorithme de décomposition par les quantités. La fonction prend en entrée les données de l'instance ainsi qu'un jeu de paramètres et renvoie en sortie la solution obtenue (à chaque iteration), le nombre d'itérations effectuées, le temps d'exécution, la valeur objective en la solution obtenue J(P*), l'ensemble des solutions calculées à chaque itération et les valeurs du multiplicateur à chaque itération.

	- DecompositionPrediction.m: Implémentation de l'algorithme de décomposition par prédiction. La fonction prend en entrée les données de l'instance ainsi qu'un jeu de paramètres et renvoie en sortie la solution obtenue (à chaque iteration), le nombre d'itérations effectuées, le temps d'exécution, la valeur objective en la solution obtenue J(P*), l'ensemble des solutions calculées à chaque itération et les valeurs du multiplicateur à chaque itération.
	
	- ResolutionExact.m: Fonction permettant de résoudre une instance données en paramètre à l'aide du solveur de MATLAB. Elle prend en entrée les données de l'instance à traiter et renvoie la solution, la valeur objective et la valeur du multiplicateur. Cette fonction nous sert de référence.

	- TestKKT.m: Fonction permettant de tester les différentes conditions KKT dans le cadre de ce problème (la fonction est générique, aussi elle prend en paramètre des multiplicateurs lambda en plus des mu ; qui sont donc à 0 dans l'exécution).



/Figures : répertoire des figures présentées dans le rapport. 