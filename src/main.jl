include("model.jl")

path = "../data/n250_p10_L20_m25_Rcouv30_Rcom30.txt"

println("=== Solving with MTZ")
solveMTZ(path, texOutput = true, isRelaxation=true)
#=
println("\n=== Solving with subtours at the root then MTZ")
solveStMTZ(path)

println("\n=== Generalized subtours at the root then MTZ")
solveStMTZG(path)

println("\n=== Subtours in callback")
solveCB(path)

println("\n=== Generalized subtours incallback")
solveCBG(path)

=#
