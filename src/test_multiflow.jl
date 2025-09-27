include("model.jl")

instancePath = "../data/expeCAID1/set7/n100_gridLength600_rCom1500.txt"
resultsMF = solve(instancePath, p = 5, time_limit = 3000, connexity = :MultiFlow)
resultsMTZ = solve(instancePath, p = 5, time_limit = 3000, connexity = :MTZ)

println("\tMTZ\tMF")
println("Obj:\t", resultsMTZ["objective"], "\t", resultsMF["objective"])
println("LB:\t", resultsMTZ["lowerBound"], "\t", resultsMF["lowerBound"])
println("Time:\t", round(Int, resultsMTZ["resolutionTime"]), "s\t", round(Int, resultsMF["resolutionTime"]), "s")
