using AutoExpe

include("../model.jl")

function aeSolveMTZCPICmax50_itmax10(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    
    resParam.cpIneqFamilies = [getICIFamily(10)]
    resParam.cpMaxCutsCount = 50

    return aeSolve(resParam)
end
function aeSolveMTZBCvide(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.threads=1
    resParam.bcIneqFamilies = []
    return aeSolve(resParam)
end

function aeSolveMTZBCWeakerICI(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.threads=1
    resParam.bcIneqFamilies = [getWeakerICIFamily(30)]
    return aeSolve(resParam)
end

function aeSolveMF3hhCPIC(dParam::Dict{String, Any})
    resParam = ResolutionParam("MultiFlow3hh", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.cpIneqFamilies = [getICIFamily(10)]
    resParam.cpMaxCutsCount = 50
    return aeSolve(resParam)
end

function aeSolveMTZCPICmax50_itmax30(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.cpIneqFamilies = [getICIFamily(30)]
    resParam.threads=1
    resParam.cpMaxCutsCount = 50
    
    return aeSolve(resParam)
end
function aeSolveMTZCPvide(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.cpIneqFamilies = []
    resParam.threads=1
    # No cuts are added
    resParam.cpMaxCutsCount = 0
    
    return aeSolve(resParam)
end
function aeSolveMTZCPICmax100_itmax10(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.cpIneqFamilies = [getICIFamily(10)]
    resParam.cpMaxCutsCount = 100

    return aeSolve(resParam)
end

function aeSolveMTZCPICmax300_itmax10(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.cpIneqFamilies = [getICIFamily(10)]
    resParam.cpMaxCutsCount = 300

    return aeSolve(resParam)
end


function aeSolveMTZ(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    return aeSolve(resParam)
end

function aeSolveMFSourour(dParam::Dict{String, Any})
    resParam = ResolutionParam("MultiFlowSourour", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    return aeSolve(resParam)
end

function aeSolveFlow(dParam::Dict{String, Any})
    resParam = ResolutionParam("Flow", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    return aeSolve(resParam)
end

function aeSolveMTZCIIneq(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    #resParam.threads=1
    resParam.ciIneq = true
    return aeSolve(resParam)
end#not only 1 threads

function aeSolveMTZ3(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ3", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    return aeSolve(resParam)
end

function aeSolveMTZC2(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZC2", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    return aeSolve(resParam)
end

function aeSolveMTZC2CIIneq(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZC2", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.ciIneq = true
    return aeSolve(resParam)
end

function aeSolveMTZMaxCliques(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.addMaxCliques = true
    return aeSolve(resParam)
end

function aeSolveMTZTightenRelax(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.tightenRelaxation = true
    return aeSolve(resParam)
end

function aeSolveMTZPP(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.preprocessing = true
    return aeSolve(resParam)
end

function aeSolveMTZ1Thread(dParam::Dict{String, Any})
    resParam = ResolutionParam("MTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.threads = 1
    return aeSolve(resParam)
end

function aeSolveMF(dParam::Dict{String, Any})
    resParam = ResolutionParam("MultiFlow", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    return aeSolve(resParam)
end

function aeSolveMF3h(dParam::Dict{String, Any})
    resParam = ResolutionParam("MultiFlow3h", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    return aeSolve(resParam)
end

function aeSolveMF3hh(dParam::Dict{String, Any})
    resParam = ResolutionParam("MultiFlow3hh", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    return aeSolve(resParam)
end

function aeSolveMF3hhCIIneq(dParam::Dict{String, Any})
    resParam = ResolutionParam("MultiFlow3hh", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.ciIneq = true
    return aeSolve(resParam)
end


function aeSolveMF3hhMC(dParam::Dict{String, Any})
    resParam = ResolutionParam("MultiFlow3hh", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.addMaxCliques = true
    return aeSolve(resParam)
end

function aeSolveMF3(dParam::Dict{String, Any})
    resParam = ResolutionParam("MultiFlow3", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    return aeSolve(resParam)
end

function aeSolveCB(dParam::Dict{String, Any})
    resParam = ResolutionParam("SubtoursCallback", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    return aeSolve(resParam)
end

function aeSolveCBPP(dParam::Dict{String, Any})
    resParam = ResolutionParam("SubtoursCallback", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.preprocessing = true
    return aeSolve(resParam)
end 

function aeSolveRootCBThenMTZ(dParam::Dict{String, Any})
    resParam = ResolutionParam("RelaxSubtoursThenMTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])

    include(dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.maxCuts = 2 * m
    return aeSolve(resParam)
end 

function aeSolveCBIF(dParam::Dict{String, Any})
    resParam = ResolutionParam("SubtoursCallback", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.maxCutsCount = 10
    return aeSolve(resParam)
end

function aeSolveQuad(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.warmStart = true
    return aeSolve(resParam)
end

function aeSolveQuadMTZ(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    return aeSolve(resParam)
end

function aeSolveQuadMF(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMF3hh", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.warmStart = true
    return aeSolve(resParam)
end

function aeSolveQuadMFCI(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMF3hh", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.cpIneqFamilies = [getICIFamily(10)]
    resParam.cpMaxCutsCount = 100
    return aeSolve(resParam)
end

function aeSolveQuadMFCIMC(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMF3hh", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.cpIneqFamilies = [getICIFamily(10), getCliqueIFamily(10)]
    resParam.cpMaxCutsCount = 100
    return aeSolve(resParam)
end

function aeSolveQuadBC_MFCIMC(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMF3hh", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.bcIneqFamilies = [getWeakerICIFamily(2), getCliqueIFamily(10)]
    return aeSolve(resParam)
end

function aeSolveQuadBCCP_MFCIMC(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMF3hh", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.bcIneqFamilies = [getICIFamily(10), getCliqueIFamily(10)]
    resParam.cpIneqFamilies = [getICIFamily(10), getCliqueIFamily(10)]
    resParam.cpMaxCutsCount = 100
    return aeSolve(resParam)
end

function aeSolveQuadCB(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadCliqueCallback", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.warmStart = true
    resParam.maxCutsCount = 10
    return aeSolve(resParam)
end

function aeSolveQuadTightenRelax(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.warmStart = true
    resParam.tightenRelaxation = true
    return aeSolve(resParam)
end

function aeSolveQuadMTZAddMaxCliques(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.cpIneqFamilies = [getCliqueIFamily(10)]
    resParam.cpMaxCutsCount = 100
    return aeSolve(resParam)
end

function aeSolveQuadMTZCI(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.cpIneqFamilies = [getICIFamily(10)]
    resParam.cpMaxCutsCount = 100
    return aeSolve(resParam)
end

function aeSolveQuadMTZMCCI(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.cpIneqFamilies = [getICIFamily(10), getCliqueIFamily(10)]
    resParam.cpMaxCutsCount = 100
    return aeSolve(resParam)
end

function aeSolveQuadMFAddMaxCliques(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMF3hh", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.cpIneqFamilies = [getCliqueIFamily(10)]
    resParam.cpMaxCutsCount = 100
    return aeSolve(resParam)
end

function aeSolveQuadTightenRelax1(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.warmStart = true
    resParam.tightenRelaxation = true
    resParam.redundancyRatio = 1.0
    return aeSolve(resParam)
end

function aeSolveQuadTightenRelax075(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.warmStart = true
    resParam.tightenRelaxation = true
    resParam.redundancyRatio = 0.75
    return aeSolve(resParam)
end

function aeSolveQuadTightenRelax05(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.warmStart = true
    resParam.tightenRelaxation = true
    resParam.redundancyRatio = 0.5
    return aeSolve(resParam)
end

function aeSolveQuadTightenRelax025(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMTZ", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.warmStart = true
    resParam.tightenRelaxation = true
    resParam.redundancyRatio = 0.25
    return aeSolve(resParam)
end

function aeSolveQuadWS(dParam::Dict{String, Any})
    resParam = ResolutionParam("QuadMTZWS", dParam["instancePath"], p = dParam["p"], time_limit = dParam["time_limit"])
    resParam.warmStart = true
    return aeSolve(resParam)
end

"""
Function that calls the solve() method from model.jl in a format suitable for the autoexpe package
(i.e., take a dictionary as an input and returns a dictionary)
"""
function aeSolve(resParam::ResolutionParam)
    
    aeResults = Dict{String, Any}() 
    
   silent = resParam.verbose == 0
    if !occursin("Quad", string(resParam.connexity))

        if !resParam.isRelaxation
            
            if !silent
	        println("=== Solve the integer linear model")
            end 
            aeResults = solve(resParam)
        end 
        
        if !silent
            println("=== Solve the continuous linear model")
        end

        tempIsRelaxation = resParam.isRelaxation
        resParam.isRelaxation = true
        relaxResults = solve(resParam)
        resParam.isRelaxation = tempIsRelaxation

        if haskey(relaxResults, "objective")
            aeResults["rootRelaxation"] = relaxResults["objective"]
        end
    else    

        if !resParam.isRelaxation

            if !silent  
                println("=== Solve the integer quadratic model")
            end 
            aeResults = solveQuadratic(resParam)
        end 

        

        println("Warning: computation ofthe relaxation in quadratic models removed")
        #=
        if !silent
            println("=== Solve the continuous quadratic model")
        end
        
        tempIsRelaxation = resParam.isRelaxation
        resParam.isRelaxation = true
        relaxResults = solveQuadratic(resParam)
        resParam.isRelaxation = tempIsRelaxation
        
        if haskey(relaxResults, "objective")
            aeResults["rootRelaxation"] = relaxResults["objective"]
        end 
=#
    end 
        
    return aeResults
end 

function aeSolveGreedy(resParam::ResolutionParam; isDeterministic::Bool=false)
    return aeSolveMetaGreedy(resParam, solveGreedy, isDeterministic = isDeterministic,p=resParam.p)
end

function aeSolveGreedyRelax(resParam::ResolutionParam; isDeterministic::Bool=false)
    return aeSolveMetaGreedy(resParam, solveGreedyRelax, isDeterministic = isDeterministic,p=resParam.p)
end

function aeSolveGreedyDeterministic(resParam::ResolutionParam)
    return aeSolveGreedy(resParam, isDeterministic = true,p=resParam.p)
end

function aeSolveGreedyRelaxDeterministic(resParam::Dict{String, Any})
    return aeSolveGreedyRelax(resParam, isDeterministic = true,p=resParam.p)
end 

"""
Use a greedy algorithm to solve an instance of the problem.
If the deterministic version of the algorithm is considered, only solve it once.
If the stochastic version is considered, solve it until:
- the time limit is reached; or
- all clients are covered; or
- 20 iterations occurred and no client is covered; or
- 20 iterations occurred without improving the solution.

Input:
- param: parameters (instance path, is it stochastic, time limit, is it connected, ...)
- greedyFunction: function that will be called to solve the problem
"""
function aeSolveMetaGreedy(resParam::ResolutionParam, greedyFunction; isDeterministic::Bool=false,p::Int=5)

    instance = Instance(resParam.instancePath)
    instance.p = p
    # If the time_limit is -1, the algorithm is only applied once without any stochasticity
    isStochastic = resParam.time_limit != -1 && !isDeterministic 

    isConnected = resParam.connexity != :None

    # Best objective values obtained accross time
    bestObjectives = Vector{Int}()
    bestObjectivesTime = Vector{Float64}()
    isSiteOpened = nothing
    isClientCovered = nothing
    clientCoveredCount = 0

    iterationCount = 0

    startingTime = time()
    
    # Stop if the best solution is not improved when this number of iterations is reached
    # (updated each time an improvment is made)
    maximalIterationsWithoutImprovement = 20
    currentResults = nothing 

    # Solve using the heuristic while:
    # - all clients ar not covered yet; and
    # - the time is not over; and
    # - the maximal number of iterations without improvement is reached; and
    # - 100 iteration have not occured without covering any client; and
    # - there is no solution or the resolution is stochastic (if it is stochastic, there is only one resolution)
    while clientCoveredCount < instance.n && (resParam.time_limit == -1 || time() - startingTime < resParam.time_limit) && iterationCount < maximalIterationsWithoutImprovement && (iterationCount < 20 || bestObjectives[end] > 0) && (length(bestObjectivesTime) == 0 || isStochastic)

        iterationCount += 1
        elapsedTime = time() - startingTime
        remainingTime = resParam.time_limit - elapsedTime

        if resParam.time_limit == -1
            remainingTime = -1
        end 
        println("Warning: client height ignored (not yet implemented in quadratic formulations)")
        if iterationCount == 1
            currentResults = greedyFunction(resParam.instancePath, isConnected = isConnected, isStochastic = false, p = instance.p, time_limit = round(Int, remainingTime))
        else
            currentResults = greedyFunction(resParam.instancePath, isConnected = isConnected, isStochastic = isStochastic, p = instance.p, time_limit = round(Int, remainingTime))
        end
        
        # If the solution obtained is better than the currently best known
        if length(bestObjectivesTime) == 0 || currentResults["objective"] > bestObjectives[end]
            clientCoveredCount = currentResults["objective"]
            elapsedTime = time() - startingTime
            
            push!(bestObjectives, clientCoveredCount)
            push!(bestObjectivesTime, elapsedTime)
            isSiteOpened = currentResults["isSiteOpened"]
            isClientCovered = currentResults["isClientCovered"]

            # Get the average time of an iteration. The algorithms will stop if this duration * 100 is reached without any improvement
            maximalIterationsWithoutImprovement = iterationCount + 20
        end
    end

    results = Dict{String, Any}()
    results["resolutionTime"] = time() - startingTime
    results["iterationCount"] = iterationCount
    results["bestObjectives"] = bestObjectives
    results["bestObjectivesTime"] = bestObjectivesTime
    results["isConnected"] = isConnected
    results["isClientCovered"] = isClientCovered
    results["isSiteOpened"] = isSiteOpened
    results["objective"] = bestObjectives[end]
    results["resolutionTimeToBest"] = bestObjectivesTime[end]
    results["instancePath"] = resParam.instancePath
    results["p"] = instance.p
    results["n"] = instance.n
    results["m"] = instance.m
    return results
end
using CSV
using DataFrames
using Statistics: argmax
"""
Batch processing function
- instance_paths: Vector{String}, each element is a data file path
- p: number of HAPS
- time_limit: time limit for each model (seconds)
- csv_path: output csv file path
"""
function choosewindowsize(instance_paths::Vector{String}; p::Int=5, time_limit::Int=-1, csv_path::String="D:/locHAPS_stage/results/resultset1.csv")
    results = DataFrame(
        instance = String[],
        windowSize = Int[],
        objective = Float64[],
        resolutionTime = Float64[],
        gap = Float64[],
        n = Int[],
        m = Int[]
    )

    for path in instance_paths
        param = ResolutionParam("QuadMTZ", path, p=p, time_limit=time_limit)
        param.warmStart = true
        param.useGridInitialization = false
        param.threads=20
        for ws in [50, 100, 150, 200, 250, 300, 400, 500, 600, 800, 1000, 1500]
            res = solveQuadratic(param; windowSize=ws)
            obj = get(res, "objective", NaN)
            time = get(res, "resolutionTime", NaN)
            lb = get(res, "lowerBound", obj)
            gap = obj == 0 ? 0.0 : abs(obj - lb) / abs(obj)
            n = get(res, "n", NaN)
            m = get(res, "m", NaN)
            push!(results, (path, ws, obj, time, gap, n, m))
        end
    end

    CSV.write(csv_path, results;delim=';')
end
function batch_solve_and_export(instance_paths::Vector{String}; p::Int=5, time_limit::Int=600, csv_path::String="D:/locHAPS_stage/results/resultset2.csv")
    #choosewindowsize(instance_paths; p=p, time_limit=time_limit)
    all_ws_df = CSV.read("D:/locHAPS_stage/results/resultset1.csv", DataFrame)
    Results = DataFrame(
        instance = String[],
        method = String[],
        objective = Float64[],
        resolutionTime = Float64[],
        gap = Float64[],
        n = Int[],
        m = Int[]
    )

    for (i, path) in enumerate(instance_paths)
        rng = MersenneTwister(1234 + i)
        row = filter(r -> r.instance == path, all_ws_df)
        if nrow(row) == 0
            println("No best window size found for instance: $path, skipping solveQuadratic with best ws.")
            continue
        end
        idx = argmax(row.objective)
        best_ws = row.windowSize[idx]
        # 1. aeSolve
        println("Computing aeSolve")
        param1=ResolutionParam("MTZ", path, p=p, time_limit=time_limit)
        param1.threads=20
        res1 = aeSolve(param1)
        obj1 = get(res1, "objective", NaN)
        time1 = get(res1, "resolutionTime", NaN)
        lb1 = get(res1, "lowerBound", obj1)
        gap1 = obj1 == 0 ? 0.0 : abs(obj1 - lb1) / abs(obj1)
        n1 = get(res1, "n", NaN)
        m1 = get(res1, "m", NaN)
        push!(Results, (path, "aeSolve", obj1, time1, gap1, n1, m1))

        # 2. aeSolveGreedy
        println("Computing aeSolveGreedy")
        res2 = aeSolveGreedy(ResolutionParam("MTZ", path, p=p, time_limit=time_limit); isDeterministic=false)
        obj2 = get(res2, "objective", NaN)
        time2 = get(res2, "resolutionTime", NaN)
        n2 = get(res2, "n", NaN)
        m2 = get(res2, "m", NaN)
        push!(Results, (path, "aeSolveGreedy", obj2, time2, NaN, n2, m2))

        # 3. solveQuadratic with warmStart=true, useGridInitialization=false
        println("Computing solveQuadratic1")
        param3 = ResolutionParam("QuadMTZ", path, p=p, time_limit=time_limit)
        param3.warmStart = true
        param3.useGridInitialization = false
        param3.threads=20
        res3 = solveQuadratic(param3;windowSize=best_ws)
        obj3 = get(res3, "objective", NaN)
        time3 = get(res3, "resolutionTime", NaN)
        lb3 = get(res3, "lowerBound", obj3)
        gap3 = obj3 == 0 ? 0.0 : abs(obj3 - lb3) / abs(obj3)
        n3 = get(res3, "n", NaN)
        m3 = get(res3, "m", NaN)
        push!(Results, (path, "solveQuadratic_warmStart", obj3, time3, gap3, n3, m3))

        # 4. solveQuadratic with warmStart=true, useGridInitialization=true
        println("Computing solveQuadratic2")
        param4 = ResolutionParam("QuadMTZ", path, p=p, time_limit=time_limit)
        param4.warmStart = true
        param4.useGridInitialization = true
        param4.threads=20
        res4 = solveQuadratic(param4;windowSize=best_ws)
        obj4 = get(res4, "objective", NaN)
        time4 = get(res4, "resolutionTime", NaN)
        lb4 = get(res4, "lowerBound", obj4)
        gap4 = obj4 == 0 ? 0.0 : abs(obj4 - lb4) / abs(obj4)
        n4 = get(res4, "n", NaN)
        m4 = get(res4, "m", NaN)
        push!(Results, (path, "solveQuadratic_warmStart_grille", obj4, time4, gap4, n4, m4))
    end

    CSV.write(csv_path, Results; delim=';')
    println("Results written to $csv_path")
end