# Example of commands to solve a problem with a formulation
jiang@morey:/home/uma/jiang$ useuma julia
jiang@morey:/home/uma/jiang$ julia
instance_paths=[
    "../data/expeCAID1/set5/n100_gridLength200_rCom1500.txt",
    "../data/expeCAID1/set5/n100_gridLength300_rCom1500.txt",
    "../data/expeCAID1/set5/n100_gridLength375_rCom1500.txt",
    "../data/expeCAID1/set5/n100_gridLength400_rCom1500.txt",
    "../data/expeCAID1/set5/n100_gridLength500_rCom1500.txt",
    "../data/expeCAID1/set5/n100_gridLength600_rCom1500.txt",
    "../data/expeCAID1/set5/n100_gridLength750_rCom1500.txt",
    "../data/expeCAID1/set5/n100_gridLength1000_rCom1500.txt",
    "../data/expeCAID1/set5/n100_gridLength1200_rCom1500.txt",
    "../data/expeCAID1/set5/n1000_gridLength200_rCom1500.txt",
    "../data/expeCAID1/set5/n1000_gridLength300_rCom1500.txt",
    "../data/expeCAID1/set5/n1000_gridLength375_rCom1500.txt",
    "../data/expeCAID1/set5/n1000_gridLength400_rCom1500.txt",
    "../data/expeCAID1/set5/n1000_gridLength500_rCom1500.txt",
    "../data/expeCAID1/set5/n1000_gridLength600_rCom1500.txt",
    "../data/expeCAID1/set5/n1000_gridLength750_rCom1500.txt",
    "../data/expeCAID1/set5/n1000_gridLength1000_rCom1500.txt",
    "../data/expeCAID1/set5/n1000_gridLength1200_rCom1500.txt"
]
cd("D:/locHAPS_stage/src")
include("autoExpe/autoExpe.jl")
instancePath = "../data/expeCAID1/set5/n100_gridLength500_rCom1500.txt"
# Type de formulation considérée
connexity = "MTZ" # Autres valeurs possibles "None", "MTZ", "MultiFlow", "MultiFlow3h", "MultiFlow3hh", "SubtoursCallback", "RelaxSubtoursThenMTZ", "QuadCliqueCallback", "QuadMTZ", "QuadMTZWS", "MTZC2", "MTZ3", "QuadMF3hh", "Flow", "MultiFlowSourour"

resParam = ResolutionParam(connexity, instancePath, p=5, time_limit=300)


# Autres paramètres fixables
resParam.displaySolution = false # true pour afficher la valeur des variables non nulles
resParam.tightenRelaxation = false # ajoute itérativement des inégalités de clique pour couper les relaxations successives à la racine
resParam.addMaxCliques = false # ajoute des inégalités de clique initialement (avant calcul de la relaxation)
resParam.preprocessing = false # retire des positions dominées
resParam.warmStart = true # utilise la solution d'un algorithme glouton en entrée des formulations
resParam.maxCutsCount = 10 
resParam.threads = -1 # limite le nombre de threads (-1 si pas de limitation)
resParam.redundancyRatio = 0.5 # Redondance des sommets dans les inégalités de clique (utilisé uniquement si ou utilise tightenRelaxation)
resParam.isRelaxation = false
resParam.useGridInitialization = false
resParam.addEigenvalue= false # Ajoute des contraintes d'eigenvalue pour la relaxation
# Résolution exacte
results = aeSolve(resParam)
texDocument(results)

# Résolution gloutonne
aeSolveGreedy(resParam, isDeterministic = false) # Basé sur le nombre de clients couverts
aeSolveGreedyRelax(resParam, isDeterministic = false) # Basé sur la relaxation du problème

texDocument(Instance(resParam.instancePath))

