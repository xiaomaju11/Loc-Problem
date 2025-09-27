# A ouvrir dans ./src/
using JSON
include("struct/instance.jl")
include("texOutput.jl")

function main()
    path = "/home/zach/Test/res/resultsCAID4_max120/res"
    outputFile = "drawing.tex"

    @show outputFile
    f = open(outputFile, "w") do f 

        write(f, "\\documentclass{article} \n \\usepackage{tikz} \n \\begin{document} \n ") 

        for file in readdir(path)

            if occursin(".json", file)
                stringdata=join(readlines(path * "/" * file))                                                                        
                data = JSON.parse(stringdata)

                for t in data
                    if haskey(t, "resolutionMethodName") && t["resolutionMethodName"] == "aeSolveMTZ"
                        instancePath = replace(t["instancePath"], "../.." => "..")
                        @show instancePath
                        instance = Instance(instancePath)
                        result = simpleTexPicture(instance, clientCoverage=Vector{Float64}(t["isClientCovered"]), siteOpening=Vector{Float64}(t["isSiteOpened"]))
                        write(f, result)
                        write(f, "\\begin{center}" * replace(file, "_" => "\\_") * "\\end{center}")
                    end 
                end 
            end
        end 
        write(f, "\\end{document} \n")
    end 
end

function mainSPGenerateInstances()
    
    outputFile = "/home/zach/Test/schemasLocHAPS/drawing.tex"

    coveredUnits, coveredUnitsCount, isSiteOpened, quadSitesXLocation, quadSitesYLocation, instancePaths = generateSmallExample()

    f = open(outputFile, "w") do f 

        println("!!!!!!!")
        write(f, "\\documentclass{article} \n \\usepackage{tikz} \n\\usetikzlibrary{backgrounds}\n \\pagestyle{empty} \n \\begin{document} \n ")
        
        for solutionId in 1:length(coveredUnits)

            currentCoveredUnitsCount = coveredUnitsCount[solutionId]
            currentCoveredUnits = coveredUnits[solutionId]
            
            result = nothing

            instance = Instance(instancePaths[solutionId])
            
            # If the solution corresponds to the solution of a discrete problem
            if solutionId < length(coveredUnits)
                result = simpleTexPicture(instance, clientCoverage=Vector{Float64}(currentCoveredUnits), siteOpening=Vector{Float64}(isSiteOpened[solutionId]))  
            else
                result = simpleTexPicture(instance, clientCoverage=Vector{Float64}(currentCoveredUnits), xHAPS=quadSitesXLocation, yHAPS=quadSitesYLocation)  
            end

            result = replace(result, "begin{tikzpicture}" => "begin{tikzpicture}[show background rectangle,inner frame sep=10mm]")
            result = replace(result, "/100" => "/2")
            result = replace(result, "end{center}" => "end{center}\\newpage\n\n")
            write(f, result)
        end 
        write(f, "\\end{document} \n") 
    end 
end

"""
Generate an instance with few units but in which the objective function improves with the granularity of the grid.
The objective function should also further improve when using the quadratic formulation
"""
function generateSmallExample()

    gridLengths = [25, 20, 10]
    rCom = 25
    rCouv = 25
    L = 20
    fieldLength = 100
    saveClientSiteDistances = false
    n = 10
    
    param = Dict{String, Any}() 
    param["isRelaxation"]       = false
    param["addMaxCliques"]      = false
    param["tightenRelaxation"]  = false
    param["time_limit"]         = -1
    param["warmStart"]          = true
    param["connexity"]          = "MTZ"
    param["instancePath"]       = "instanceWithGap.txt"
    param["p"]                  = 4
    param["verbose"]            = 0

    coveredUnits = nothing
    coveredUnitsCount = nothing
    isSiteOpened = nothing
    quadSitesXLocation = nothing
    quadSitesYLocation = nothing
    instancesPath = nothing

    outputPath = "../../data/"
    
    isValid = false
    while !isValid

        isValid = true

        gridLengthId = 1

        coveredUnitsCount = Vector{Int}()
        coveredUnits = Vector{Vector{Int}}()
        isSiteOpened = Vector{Vector{Int}}()
        instancesPath = Vector{String}()
            
        generateGridInstances(n, fieldLength, gridLengths, rCom, rCouv, L, 1, outputPath, saveClientSiteDistances = saveClientSiteDistances) 

        # For each grid length and while improving the granularity improves the objective value
        while isValid && gridLengthId <= length(gridLengths)

            gridLength = gridLengths[gridLengthId]
            println("++++ Resolution for grid length ", gridLength)
            param["connexity"] = "MTZ" 

            outputFile = outputPath * "n" * string(n) * "_gridLength" * string(gridLength) * "_rCom" * string(rCom) * ".txt"
            param["instancePath"]       = outputFile
            push!(instancesPath, outputFile)

            results = aeSolve(param)
            @show gridLength, round(Int, results["objective"])
            
            if length(coveredUnits) == 0 || coveredUnitsCount[end] < results["objective"]
                push!(coveredUnitsCount, round(Int, results["objective"]))
                push!(coveredUnits, round.(Int, results["isClientCovered"]))
                push!(isSiteOpened, round.(Int, results["isSiteOpened"]))
                gridLengthId += 1
            else
                isValid = false
                gridLengthId = 1
                coveredUnitsCount = Vector{Int}()
                coveredUnits = Vector{Vector{Int}}()
                isSiteOpened = Vector{Vector{Int}}()
                instancesPath = Vector{String}() 
            end
        end

        if isValid
            
            param["connexity"] = "QuadMTZ"
            param["addMaxCliques"] = true
            resultsQuad = aeSolve(param)
            
            push!(instancesPath, instancesPath[end])

            @show "quad", round(Int, resultsQuad["objective"])
            if coveredUnitsCount[end] < resultsQuad["objective"]
                push!(coveredUnitsCount, round(Int, resultsQuad["objective"]))
                push!(coveredUnits, round.(Int, resultsQuad["isClientCovered"]))
                quadSitesXLocation = resultsQuad["xHAPS"]
                quadSitesYLocation = resultsQuad["yHAPS"]
            else
                isValid = false
                gridLengthId = 1
                coveredUnitsCount = Vector{Int}()
                coveredUnits = Vector{Vector{Int}}()
                isSiteOpened = Vector{Vector{Int}}()
                instancesPath = Vector{String}() 
            end
        end
    end

    return coveredUnits, coveredUnitsCount, isSiteOpened, quadSitesXLocation, quadSitesYLocation, instancesPath
end 

#=
instance1500 = Instance("../data/expeCAID1/set5/n100_gridLength1500_rCom1500.txt")
isSiteOpened1500 = [ -0.0, 1.0, 1.0, 0.0, -0.0, 0.0, 0.0, 1.0, 1.0, -0.0, -0.0, 1.0, 0.0, 1.0, -0.0, -0.0, 1.0, 1.0, 1.0, 0.0, -0.0, 1.0, 0.0, 0.0, -0.0]
    
isClientCovered1500 = [ 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, -0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

simpleTexPicture(instance1500, clientCoverage=isClientCovered1500, siteOpening=isSiteOpened1500, path="./1500.tex")


instance1000 = Instance("../data/expeCAID1/set5/n100_gridLength1000_rCom1500.txt")
isSiteOpened1000 = [ 0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0, -0.0, 0.0, 1.0, -0.0, 1.0, 0.0, 0.0, -0.0, -0.0, -0.0, 1.0, 0.0, -0.0, -0.0, -0.0, 1.0, 1.0, 0.0, -0.0, 0.0, -0.0, -0.0, -0.0, -0.0, 1.0, 0.0, 1.0, 0.0, -0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, -0.0, 1.0, -0.0, 0.0, 0.0, -0.0, -0.0]

isClientCovered1000 = [ 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -0.0, -0.0, 1.0, 1.0, 0.0, 1.0, -0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, -0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]


simpleTexPicture(instance1000, clientCoverage=isClientCovered1000, siteOpening=isSiteOpened1000, path="./1000.tex")


instance500 = Instance("../data/expeCAID1/set5/n100_gridLength500_rCom1500.txt")
isSiteOpened500 = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	

isClientCovered500 = [ 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

simpleTexPicture(instance500, clientCoverage=isClientCovered500, siteOpening=isSiteOpened500, path="./500.tex")
=#
