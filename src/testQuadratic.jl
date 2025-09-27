include("model.jl")

timelimit = 3000
instancePath = "../data/expeCAID1/set7/n100_gridLength600_rCom1500.txt"

resultsQuadratic = solveQuadratic(instancePath, p = 5, time_limit = timelimit, useWarmStart = false)
resultsQuadraticWM = solveQuadratic(instancePath, p = 5, time_limit = timelimit, useWarmStart = true)
resultsMTZ = solve(instancePath, p = 5, time_limit = timelimit, connexity = :MTZ)

texDocument(resultsQuadratic, path="../results/resultsQ.tex")
texDocument(resultsQuadraticWM, path="../results/resultsQWM.tex")
texDocument(resultsMTZ, path="../results/resultsMTZ.tex")

println("\tMTZ\tQuad\tQuadWM")
println("Clients:\t", resultsMTZ["objective"], "\t", resultsQuadratic["objective"], "\t", resultsQuadraticWM["objective"])
println("LB:\t", resultsMTZ["lowerBound"], "\t", resultsQuadratic["lowerBound"], "\t", resultsQuadraticWM["lowerBound"])
println("Time:\t", round(Int, resultsMTZ["resolutionTime"]), "s\t", round(Int, resultsQuadratic["resolutionTime"]), "s\t", round(Int, resultsQuadraticWM["resolutionTime"]))
