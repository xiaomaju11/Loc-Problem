# Premier étape du stage

Le premier objectif du stage est de vérifier que le code fonctionne correctement et permet de résoudre des instances.
Nous avons défini plusieurs modélisations du problème. La plus simple est appellée MTZ. Pour résoudre une instance avec ce modèle, voici la procédure à suivre (plus de détails dans le fichier src/00_what_to_do.jl).

1. se  déplacer dans  le répertoire  `src`,
2. ouvrir une console julia dans ce répertoire
3. Exécuter les commandes suivantes

    include("autoExpe/autoExpe.jl") # Charge les fonctions julia qu'on va utiliser
    param = Dict{String, Any}([]) # Dictionnaire qui définit quelle instance on résout et comment
    param["connexity"] = "MTZ" # Indique qu'on utilise la formulation MTZ du problème pour résoudre 
    param["instancePath"] = "../data/expeCAID1/set7/n100_gridLength1000_rCom1500.txt" # Chemin du fichier contenant l'instance qu'on va résoudre
    param["p"] = 2 # Nombre de drones
    results = aeSolve(param) # Lancement de la résolution

Pour avoir une représentation graphique d'une solution, on  peut
utiliser :

    texDocument(results)
	
# Organisation des fichiers

*  `src/struct/instance.jl`   :  contient  la  définition   de  la  structure
  `Instance` qui regroupe toutes les informations d'une instance. On y
  trouve   également  la   fonction   `generateInstance`  qui   génère
  aléatoirement une instance en fonction des paramètres d'entrée ;
*  `src/model.jl` :  contient  les fonctions  permettant  de créer  et
  résoudre les modèles correspondant aux différentes modélisations ;
* `src/texOutput.jl` : contient les  fonctions permettant de créer des
  fichiers `.tex` représentant des solutions ;
* `data` : contient les instances (1 fichier par instance) ;
*  `autoExpe` :  contient les  fichiers relatifs  aux expérimentations
  automatiques avec le package [AutoExpe.jl](https://github.com/ZacharieALES/AutoExpe.jl).

# Second étape du stage

Une fois qu'on a vérifié que le code fonctionne, l'objectif sera de développer une nouvelle méthode de résolution basé sur deux résolutions successives :
1. Résoudre le problème avec un modèle linéaire (positions possibles pour les drones discrètes)
2. Résoudre le problème avec un problème quadratique (positions possibles pour les drones continues) en imposant que les drones restent proches des positions des drones dans la solution linéaire.

Pour l'étape 1, il suffit de résoudre un modèle linéaire (ex : MTZ). Cela ne demande pas d'adaptation.
Pour l'étape 2, il faut récupérer la solution de l'étape 1 et la donner en entrée d'une méthode de résolution quadratique (par exemple QuadMTZ) et adapter la fonction getQuadraticModel du fichier model.jl pour ajouter des contraintes qui forcent les drones à être proches de la solution précédente.

On pourra aussi donner en paramètre un rayon qui indiquera de combien au maximum la position des drones peuvent s'éloigner de la solution précédente.

Ensuite comparer :
- le modèle linéaire seul ;
- le modèle quadratique seul ;
- la nouvelle méthode.
(en faisant varier le rayon, le nombre de positions considérées dans le modèle linéaire, etc.)

La nouvelle méthode devrait avoir un objectif au moins aussi bon que le modèle linéaire seul.
Le modèle quadratique ne devrait pas être très bon car si on ne contraint pas la position des drones, il est extrêmement long à converger.

# Troisième étape du stage (pour aller plus loin)

Comparer les différents modèles linéaires, quadratiques, etc.
Essayer de faire fonctionner les approches basées sur des callbacks (i.e., génération dynamique d'inégalités durant la résolution).