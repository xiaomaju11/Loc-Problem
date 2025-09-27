using JuMP, Random, CPLEX, StatsBase
include("struct/instance.jl")
include("texOutput.jl")

"""
Solve the model with MTZ 

Input:
- instancePath: path of the instance file
- (optional) time_limit: resolution time limit (no limit if = -1) (default -1)
- (optional) isRelaxation: true if the linear relaxation is solved (default false)
- (optional) texOutput: true if a tex file which represents the solution is created
"""
function solveMTZ(instancePath::String; time_limit::Int=-1, isRelaxation::Bool=false, texOutput::Bool=false)
    return solve(instancePath, connexity = :MTZ, time_limit=time_limit, isRelaxation = isRelaxation)
end

"""
Solve the model by generating subtour cuts at the root and the using with MTZ 

Input:
- instancePath: path of the instance file
- time_limit: resolution time limit (no limit if = -1) (default -1)
- isRelaxation: true if the linear relaxation is solved (default false)
- (optional) texOutput: true if a tex file which represents the solution is created
"""
function solveStMTZ(instancePath::String; time_limit::Int= -1, isRelaxation::Bool=false, texOutput::Bool=false)
    return solve(instancePath, connexity = :RelaxSubtoursThenMTZ, time_limit = time_limit, isRelaxation = isRelaxation)
end

"""
Solve the model by generating subtour cuts in a callback

Input:
- instancePath: path of the instance file
- (optional) time_limit: resolution time limit (-1 if no limit) (default -1)
- (optional) isRelaxation: true if the linear relaxation is solved (default false)
- (optional) texOutput: true if a tex file which represents the solution is created
"""
function solveCB(instancePath::String; time_limit::Int= -1, isRelaxation::Bool=false, texOutput::Bool=false)
    return solve(instancePath, connexity = :SubtoursCallback, time_limit = time_limit, isRelaxation = isRelaxation)
end 

"""
General method to solve the localisation problem

Input
- instancePath: path of the instance to solve
- (optional) connexity: how the connexity is handled:
  - :None : No connexity imposed (default value)
  - :MTZ : Use MTZ modelisation
  - :RelaxSubtoursThenMTZ : Use subtours elimination on the relaxation and then add MTZ to solve as integer
  - :SubtoursCallback : Generate subtour elimination cuts in a callback
- (optional) time_limit: resolution time limit (-1 if no limit) (default -1)
- (optional) maxCuts: maximal number of subtour cuts added at the root (used only if connexity is :RelaxSubtoursThenMTZ) (default 50)
- (optional) texOutput: true if a tex file which represents the solution is created
- (optional) tightenRelaxation: true if the root relaxation is improved by the addition of clique inequalities
- (optional) redundancyRatio: if the relaxation is tightened by clique inequalities, it corresponds to the percentage in [0, 1] of units in a new clique that must not have appeared in any previous clique
  = 0 : no constraint
  = 0.5 : half the units in a clique must not have appeared in a previous clique
  = 1 : all the units must not have appeared in a previous clique
"""
function solve(instancePath::String; connexity::Symbol=:None, time_limit::Int=-1, maxCuts::Int=50, isRelaxation::Bool=false, texOutput::Bool=false, p::Int=-1, maxCutsCount::Int=0, threads::Int=-1, preprocessing::Bool=false, tightenRelaxation::Bool=false, redundancyRatio::Float64=Float64(0.0))
    
    print("Reading the instance...")
    startingReading = time()
    instance = Instance(instancePath)
    println(" done in ", round(Int, time() - startingReading), "s")

    if p != -1
        instance.p = p
    end
    
    startingTime = time()
    remainingTime = time_limit

    if !isRelaxation
        print("Computing greedy solution...")
        greedyResults = solveGreedy(instance, isConnected = connexity != :None, isStochastic = false, p = instance.p, time_limit = time_limit)

        elapsedTime = greedyResults["resolutionTime"]
        greedyIsSiteOpened = greedyResults["isSiteOpened"]
        greedyIsClientCovered = greedyResults["isClientCovered"]
        greedyObjective = greedyResults["objective"]
        println(" Found solution of value ", greedyObjective, " in ", round(Int, time() - startingTime), "s")
        remainingTime = time_limit - elapsedTime 
    end 


    if time_limit == -1
        remainingTime = -1
    end

    ppResults = nothing
    if preprocessing
        ppResults = doPreprocessing(instance)
    end 

    model = nothing
    x = nothing
    y = nothing
    z = nothing
    f = nothing

    ## Get the model
    rootResults = nothing

    creationTime = time()
    print("Creating the model...")
    # If subtour inequalities are solved only at the root
    if connexity == :RelaxSubtoursThenMTZ
        model, x, y, z, rootResults = generateRelaxedSubtours(instance,  maxCuts, round(Int, remainingTime), isRelaxation = isRelaxation, threads = threads, upperBound = round(Int, greedyObjective), ppResults=ppResults)
    else
        model, x, y, z, f, r = getModel(instance, connexity = connexity, time_limit = round(Int, remainingTime), isRelaxation = isRelaxation, maxCutsCount = maxCutsCount, threads = threads, ppResults=ppResults)
    end
    println(" done in ", round(Int, time() - creationTime), "s")

    results = Dict{String, Any}()
    
    if tightenRelaxation && !isRelaxation
        tightenTime = time()
        println("Tighten the relaxation...")
        cliques, firstRelaxation, lastRelaxation = tightenRelaxationByCliques(instance, connexity = connexity, p = instance.p, redundancyRatio=redundancyRatio, time_limit = time_limit)
        println("done in ", round(Int, time() - tightenTime), "s with ", length(cliques), " inequalities")

        @constraint(model, [clique in cliques], sum(z[u] for u in clique) <= instance.p)

        results["initialBound"] = firstRelaxation
        results["tightenedBound"] = lastRelaxation
        results["tighteningTime"] = time() - tightenTime
        results["inequalitiesCount"] = length(cliques)
        
    end 

    # Set the warm start by providing the greedy solution
    if !isRelaxation
        vars=all_variables(model)
        firstYIndex = findfirst(vars .== y[1])
        set_start_value.(vars[firstYIndex:firstYIndex + instance.m - 1], round.(Int, greedyIsSiteOpened))
    end 
    
    results["instancePath"] = instancePath
    results["p"] = instance.p
    results["n"] = instance.n
    results["m"] = instance.m
    results["maxCutsCount"] = maxCutsCount

    if ppResults != nothing
        results["dominated"] = length(ppResults["dominated"])
        results["comPriorities"] = size(ppResults["comPriorities"], 1)
        results["priorities"] = size(ppResults["priorities"], 1)
    end

    if remainingTime != -1
        remainingTime = time_limit - (time() - startingTime)
    end

    if remainingTime >= 5 || remainingTime == -1
    
        set_optimizer_attribute(model, "CPX_PARAM_TILIM", remainingTime)
	
        ## Solve it
        print("Solve the model...")
        optimize!(model)
        println("done")
        resolutionTime = time() - startingTime
    
        ## Get the results
        results["resolutionTime"] = resolutionTime
    
        if !isRelaxation
            results["nodeCount"] = JuMP.node_count(model)
        end 
    
        if rootResults != nothing
            merge!(results, rootResults)
        end 
    
        if primal_status(model) == MOI.FEASIBLE_POINT
            
            results["isFeasible"] = true
            results["objective"] = JuMP.objective_value(model)
            results["isSiteOpened"] = JuMP.value.(y)
            results["isClientCovered"] = JuMP.value.(z)
            
            if termination_status(model) == MOI.OPTIMAL
                results["isOptimal"] = true
                results["lowerBound"] = results["objective"]
            else
                results["isOptimal"] = false
                results["lowerBound"] = JuMP.objective_bound(model)
            end  
        else
            println("No feasible solution found")
            results["isFeasible"] = false
	end
    else
        results["isOptimal"] = false
        results["objective"] = greedyObjective
    end

    return results
end 


"""
Solve the quadratic formulation of the localisation problem

Input
- instancePath: path of the instance to solve
- (optional) time_limit: resolution time limit (-1 if no limit) (default -1)
- (optional) texOutput: true if a tex file which represents the solution is created
- (optional) tightenRelaxation: true if the root relaxation is improved by the addition of clique inequalities
- (optional) redundancyRatio: if the relaxation is tightened by clique inequalities, it corresponds to the percentage in [0, 1] of units in a new clique that must not have appeared in any previous clique
  = 0 : no constraint
  = 0.5 : half the units in a clique must not have appeared in a previous clique
  = 1 : all the units must not have appeared in a previous clique

"""
function solveQuadratic(instancePath::String; time_limit::Int=-1, isRelaxation::Bool=false, texOutput::Bool=false, p::Int=-1, maxCutsCount::Int=0, threads::Int=-1, useWarmStart::Bool=true, tightenRelaxation::Bool=false, connexity::Symbol=:QuadMTZ, redundancyRatio::Float64=Float64(0.0))

    print("Reading the instance...") 
    startingReading = time()
    instance = Instance(instancePath)
    println(" done in ", round(Int, time() - startingReading), "s")

    if p != -1
        instance.p = p
    end
    
    startingTime = time()

    elapsedTime = 0
    
    if useWarmStart && !isRelaxation
        print("Computing greedy solution...")
        greedyResults = solveGreedy(instance, isConnected = true, isStochastic = false, p = instance.p, time_limit = time_limit)

        elapsedTime = greedyResults["resolutionTime"]
        greedyIsSiteOpened = greedyResults["isSiteOpened"]
        greedyIsClientCovered = greedyResults["isClientCovered"]
        greedyObjective = greedyResults["objective"]
        println(" Found solution of value ", greedyObjective, " in ", round(Int, time() - startingTime), "s") 
    end 
        
    remainingTime = time_limit - elapsedTime

    if time_limit == -1
        remainingTime = -1
    end

    ## Get the model    
    print("Creating the model...")
    creationTime = time()
    model, alpha, beta, z, c, t, x = getQuadraticModel(instance, time_limit = round(Int, remainingTime), isRelaxation = isRelaxation, threads = threads)
    println(" done in ", round(Int, time() - creationTime), "s")

    results = Dict{String, Any}()
    if tightenRelaxation && !isRelaxation
        tightenTime = time()
        println("Tighten the relaxation...")
        cliques, firstRelaxation, lastRelaxation = tightenRelaxationByCliques(instance, connexity = :QuadMTZ, p = instance.p, redundancyRatio=redundancyRatio, time_limit = time_limit)
        println("done in ", round(Int, time() - tightenTime), "s with ", length(cliques), " inequalities")

        @constraint(model, [clique in cliques], sum(z[u] for u in clique) <= instance.p)

        results["initialBound"] = firstRelaxation
        results["tightenedBound"] = lastRelaxation
        results["tighteningTime"] = time() - tightenTime
        results["inequalitiesCount"] = length(cliques)
        
    end 

    ## Set the warm start by providing the greedy solution
    if useWarmStart && !isRelaxation

        # Get all the opened site ordered according to the symmetry breaking constraints
        # (i.e., smallest alpha first and if same alpha, smallest beta first)
        greedyPositions = Vector{Tuple{Float64, Float64}}([])
        
        println("Haps id in greedy solution: ") 
        for j in 1:instance.m
            if greedyIsSiteOpened[j] > 1-1^-4
                push!(greedyPositions, (instance.sitesPositions[j, 1], instance.sitesPositions[j, 2]))
                println("\t", j, " (", round(instance.sitesPositions[j, 1], digits = 2), ", ", round(instance.sitesPositions[j, 2], digits = 2), ")")
            end
        end
        maxDiffY = maximum(instance.sitesPositions[:, 2]) - minimum(instance.sitesPositions[:, 2])
        sort!(greedyPositions, by = x -> x[1] * maxDiffY + x[2])

        println("Sorted positions: ", greedyPositions)

        pGreedy = length(greedyPositions)
        
        # Set the warm start variables 
        vars=all_variables(model)        
        xMin = minimum(instance.clientPositions[:, 1]) 
        xMax = maximum(instance.clientPositions[:, 1])
        yMin = minimum(instance.clientPositions[:, 2])
        yMax = maximum(instance.clientPositions[:, 2])
        
        for hId in 1:pGreedy
            h = greedyPositions[hId]
            vAlphaId = findfirst(vars .== alpha[hId])
            vBetaId = findfirst(vars .== beta[hId])
            set_start_value(vars[vAlphaId], min(xMax, max(xMin, h[1])))
            set_start_value(vars[vBetaId], min(yMax, max(yMin, h[2])))

            println("Set HAPS ", hId, " at position ", min(xMax, max(xMin, h[1])), ",", min(yMax, max(yMin, h[2])))
        end
    end
    
    results["instancePath"] = instancePath
    results["p"] = instance.p
    results["n"] = instance.n
    results["m"] = instance.m

    if remainingTime != -1
        remainingTime = time_limit - (time() - startingTime)
    end

    if remainingTime >= 5 || remainingTime == -1
    
        set_optimizer_attribute(model, "CPX_PARAM_TILIM", remainingTime)
     	
        ## Solve it
        optimize!(model)
        resolutionTime = time() - startingTime
        
        ## Get the results
        results["resolutionTime"] = resolutionTime
    
        if !isRelaxation
            results["nodeCount"] = JuMP.node_count(model)
        end 

        if primal_status(model) == MOI.FEASIBLE_POINT

            xHAPS = JuMP.value.(alpha)
            yHAPS = JuMP.value.(beta)
            
            d = zeros(instance.p, instance.p)
            for i in 1:instance.p
                for j in 1:instance.p
                    d[i, j] = sqrt((xHAPS[i]-xHAPS[j])^2+(yHAPS[i]-yHAPS[j])^2)
                    d[j, i] = d[i, j]
                end
            end
    
            results["isFeasible"] = true
            results["objective"] = JuMP.objective_value(model)
            results["xHAPS"] = JuMP.value.(alpha)
            results["yHAPS"] = JuMP.value.(beta)
            results["isClientCovered"] = JuMP.value.(z)
            
            if termination_status(model) == MOI.OPTIMAL
                results["isOptimal"] = true
                results["lowerBound"] = results["objective"]
            else
                results["isOptimal"] = false
                results["lowerBound"] = JuMP.objective_bound(model)
            end  
        else
            println("No feasible solution found")
            results["isFeasible"] = false
	end
    end

    return results
end 

"""
Get the model associated with a formulation of the localisation problem

Input
- instance: the instance solved
- connexity: see function solve()
- isRelaxation: true if the variable are all continuous
- maxCutsCount : nombre maximal de coupes ajoutées à un noeud du branch-and-bound

Output:
- model: the model
- x, y, z: its variables
"""
function getModel(instance::Instance; connexity::Symbol=:None, time_limit::Int=-1, isRelaxation::Bool=false, maxCutsCount::Int=0, threads::Int=-1, ppResults=nothing)
    
    ### Model declaration
    model = Model(CPLEX.Optimizer)
    if isRelaxation
        set_silent(model)
    end 
    if time_limit != -1 && time_limit > 5
        set_optimizer_attribute(model, "CPX_PARAM_TILIM", time_limit)
    end

    if threads != -1
        MOI.set(model, MOI.NumberOfThreads(), threads)
    end 
        
    ### Variables
    x = nothing
    r = nothing
    q = nothing
    
    if isRelaxation
        
        # 1 iff client i is covered
        @variable(model, 0 <= z[1:instance.n] <= 1)
        # 1 iff sites j is opened
        @variable(model, 0 <= y[1:instance.m] <= 1)

        # 1 iff (j,j1) is a connection 
        if !(connexity in [:None, :MultiFlow, :MultiFlow3, :MultiFlow3h, :MultiFlow3hh])
            if connexity == :SubtoursCallback
                @variable(model, 0 <= x[j in 1:instance.m, j1 in 1:instance.m ; j1>j] <= 1)
            else
                @variable(model, 0 <= x[j in 1:instance.m, j1 in 1:instance.m ; j1!=j] <= 1)
            end
        end

        if connexity == :MultiFlow3
            @variable(model, 0 <= v[j in 1:instance.m] <= 1)
        end
        
        if connexity in [:MultiFlow3h, :MultiFlow3hh]
            @variable(model, 0 <= r[j in 1:instance.m, k in 1:instance.p] <= 1)
        end 
    else
        
        # 1 iff client i is covered
	# (can be continuous but does not necessarily give better results)
        @variable(model, z[1:instance.n], Bin)
        # 1 iff sites j is opened
        @variable(model, y[1:instance.m], Bin)
        
        if !(connexity in [:None, :MultiFlow, :MultiFlow3, :MultiFlow3h, :MultiFlow3hh] )
            # x[j, j1] = 1 iff (j,j1) is a connection

            # These variables can be continuous if the subtour formulation is considered
            if connexity == :SubtoursCallback
                @variable(model, 0 <= x[j in 1:instance.m, j1 in 1:instance.m ; j1>j] <= 1)
            else
                @variable(model, x[j in 1:instance.m, j1 in 1:instance.m ; j1!=j], Bin)
            end
        end 
    
        if connexity in [:MultiFlow3h, :MultiFlow3hh]
            @variable(model, r[j in 1:instance.m, k in 1:instance.p], Bin)
        end
        
        if connexity == :MultiFlow3
            @variable(model, v[j in 1:instance.m], Bin)
        end 
    end

    f = nothing

    if connexity == :MultiFlow

        # 1 iff the arc (j1, j2) is on the path from s to t
        @variable(model, 0 <= f[s in 1:instance.m, t in 1:instance.m, j1 in 1:instance.m, j2 in 1:instance.m; isFVariable(s, t, j1, j2, instance)] <= 1)
        
    end

    if connexity == :MultiFlow3

        # 1 iff the arc (j1, j2) is on the path from the root to position j (j1 = 0 corresponds to the root)
        @variable(model, 0 <= f[j in 1:instance.m, j1 in 0:instance.m, j2 in 1:instance.m; (j1 == 0 || instance.d_com[j1,j2] <= instance.rCom) && j1 != j2 && j1 != j] <= 1)
    end 
    
    if connexity == :MultiFlow3h

        # 1 iff the arc (j1, j2) is on the path from the root to HAPS number k (j1 = 0 corresponds to the root)
        @variable(model, 0 <= f[k in 1:instance.p, j1 in 0:instance.m, j2 in 1:instance.m;  j1 == 0 || instance.d_com[j1,j2] <= instance.rCom && j1 != j2] <= 1)
    end 
    
    if connexity == :MultiFlow3hh

        # 1 iff the arc (j1, j2) is on the path from the root to HAPS number k (j1 = 0 corresponds to the root)
        @variable(model, 0 <= f[k in 2:instance.p, k1 in 1:instance.p, k2 in 1:instance.p; k1 != k2 && k1 != k] <= 1)

        @variable(model, 0 <= q[k1 in 1:instance.p, k2 in 1:instance.p; k1 != k2] <= 1)
    end 

    if ppResults != nothing
        @constraint(model, [j in ppResults["dominated"]], y[j] == 0)

        for line in 1:size(ppResults["priorities"], 1)
            j1 = ppResults["priorities"][line, 1]
            j2 = ppResults["priorities"][line, 2]
            @constraint(model, y[j2] <= y[j1])
        end

        for dId in 1:size(ppResults["comPriorities"], 1)

            # If site j2 has an HAPS, then at least one HAPS in unreachable site also has one
            d = ppResults["comPriorities"][dId]
            j1 = d["site1"]
            j2 = d["site2"]
            unreachableSites = d["unreachableSites"]

            @constraint(model, y[j2] <= sum(y[j] for j in unreachableSites if j2 != j))
        end
    end  

    if connexity == :MTZ
        # Index of site j in the connected arborescence
        @variable(model, t[j in 1:instance.m] >= 0)
    end 

    ### Objective
    @objective(model, Max, sum(z[i]  for i in 1:instance.n))

    ### Constraints
    # i cannot be covered if none of the sites within reach are opened
    @constraint(model, [i in 1:instance.n], z[i] <=  sum(y[j] for j in 1:instance.m if (instance.d[i,j]<=instance.rCouv))  )

    @constraint(model,  sum(y[j] for j in 1:instance.m) <= instance.p)

    if connexity  != :None

        # Far away sites cannot be directly connected to each other
        # Remark: already imposed in the multiflow formulation in the definition of the variables
        if !(connexity in [:MultiFlow, :MultiFlow3, :MultiFlow3h, :MultiFlow3hh])
            @constraint(model,  [j in 1:instance.m, j1 in 1:instance.m ; xExist(j, j1,  connexity) && instance.d_com[j,j1] > instance.rCom] , x[j,j1]== 0)
        end 

        if connexity == :SubtoursCallback
            
            # If j is closed, it has no connection
            @constraint(model,  [j in 1:instance.m, j1 in 1:instance.m ; j!=j1] , x[min(j,j1), max(j, j1)]<= y[j])

            # Number of connections = Number of opened sites - 1
            @constraint(model,  sum(x[j,j1] for j in 1:instance.m for j1 in 1:instance.m if xExist(j, j1, connexity)) == sum(y[j] for j in 1:instance.m) -1)
            
        elseif connexity == :MultiFlow

            # If a site is closed, it has no connection
            @constraint(model, [s in 1:instance.m, t in 1:instance.m, j1 in 1:instance.m; s < t], sum(f[s, t, j1, j2] for j2 in 1:instance.m if isFVariable(s, t, j1, j2, instance)) <= y[j1])
            
            @constraint(model, [s in 1:instance.m, t in 1:instance.m, j2 in 1:instance.m; s < t], sum(f[s, t, j1, j2] for j1 in 1:instance.m if isFVariable(s, t, j1, j2, instance)) <= y[j2])
            
            # Flow conservation: Intermediate nodes
            @constraint(model, [s in 1:instance.m, t in 1:instance.m, j1 in 1:instance.m; s < t && j1 != s && j1 != t], sum(f[s, t, j1, j2] for j2 in 1:instance.m if j2 != s && j1 != j2 && instance.d_com[j1,j2] <= instance.rCom) == sum(f[s, t, j2, j1] for j2 in 1:instance.m if j2 != t && j1 != j2 && instance.d_com[j1,j2] <= instance.rCom))
            
            # Starting node
            @constraint(model, [s in 1:instance.m, t in 1:instance.m; s < t], sum(f[s, t, s, j1] for j1 in 1:instance.m if j1 != s &&  instance.d_com[s,j1] <= instance.rCom) >= y[s] + y[t] - 1)

        elseif connexity == :MultiFlow3

            # If a site is closed, it has no connection
            @constraint(model, [j in 1:instance.m, j1 in 1:instance.m, j2 in 1:instance.m; j1 != j && j1 != j2 && instance.d_com[j1,j2] <= instance.rCom], f[j, j1, j2] <= y[j1])
            @constraint(model, [j in 1:instance.m, j1 in 0:instance.m, j2 in 1:instance.m; j1 != j && (j1 == 0 || (j1 != j2 && instance.d_com[j1,j2] <= instance.rCom))], f[j, j1, j2] <= y[j2])
            
            ## Flow conservation constraints
            # Intermediate nodes
            @constraint(model, [j in 1:instance.m, j1 in 1:instance.m; j1 != j], sum(f[j, j1, j2] for j2 in 1:instance.m if j1 != j2 && instance.d_com[j1,j2] <= instance.rCom) == sum(f[j, j2, j1] for j2 in 0:instance.m if j2 == 0 || (j1 != j2 && j2 != j && instance.d_com[j2,j1] <= instance.rCom)))
            
            ## Flow conservation: starting node

            # v[j] = 1 for only one position j
            @constraint(model,  sum(v[j] for j in 1:instance.m) == 1)
            
            # The flow of position j goes from the root to the site j1 for which v[j1] = 1
            @constraint(model, [j in 1:instance.m, j1 in 1:instance.m], f[j, 0, j1] <= v[j1])

            # One unit of the flow of position j leaves the root iff site j is opened
            @constraint(model, [j in 1:instance.m], sum(f[j, 0, j1] for j1 in 1:instance.m) == y[j])

        elseif connexity in [:MultiFlow3h, :MultiFlow3hh]

            # HAPS number 1 is affected to exactly one site
            @constraint(model, sum(r[j, 1] for j in 1:instance.m) == 1)

            # HAPS k+1 can only be affected, if HAPS k is also affected
            @constraint(model, [k in 1:instance.p-1], sum(r[j, k] for j in 1:instance.m) >= sum(r[j, k+1] for j in 1:instance.m))

            # Site j is affected an HAPS iff it is opened
            @constraint(model, [j in 1:instance.m], sum(r[j, k] for k in 1:instance.p) == y[j])

            # If HAPS number k is affected to position j, no HAPS <k can be affected to a position >j
            # (optional symmetry breaking constraints)
            @constraint(model, [j in 1:instance.m, j2 in j+1:instance.m, k in 1:instance.p, k2 in 1:k-1], 1 - r[j, k] >= r[j2, k2])

            if connexity == :MultiFlow3h
                
                # Source node flow: the flow of an affected HAPS starts by going from the root to the position of the first HAPS
                # i.e. f[k, 0, j] = r[j, 1] * sum(r[j, k] for all k)
                @constraint(model, [k in 1:instance.p, j in 1:instance.m], f[k, 0, j] <= r[j, 1])
                @constraint(model, [k in 1:instance.p, j in 1:instance.m], f[k, 0, j] <= sum(r[j, k] for j in 1:instance.m))
                @constraint(model, [k in 1:instance.p, j in 1:instance.m], f[k, 0, j] >= r[j, 1] + sum(r[j, k] for j in 1:instance.m) - 1)

                # This constraint could be used if the number of HAPS affected was equal to p rather than <= p
                #@constraint(model, [k in 1:instance.p, j in 1:instance.m], f[k, 0, j] == r[j, 1])

                # No flow on arc (j1, j2) they are not both opened
                @constraint(model, [k in 1:instance.p, j1 in 1:instance.m, j2 in 1:instance.m; j1 != j2 && instance.d_com[j1,j2] <= instance.rCom], f[k, j1, j2] <= y[j1])
                @constraint(model, [k in 1:instance.p, j1 in 0:instance.m, j2 in 1:instance.m; j1 == 0 || (j1 != j2 && instance.d_com[j1,j2] <= instance.rCom)], f[k, j1, j2] <= y[j2])

                # Flow conservation
                # if HAPS k is not in j1 (r[j1, k]=0), classical flow conservation
                # if HAPS k is in j1 (r[j1, k]=1), one more unit enters j1 than leaves it
                @constraint(model, [k in 1:instance.p, j1 in 1:instance.m], r[j1, k] + sum(f[k, j1, j2] for j2 in 1:instance.m if j1 != j2 && instance.d_com[j1,j2] <= instance.rCom) == sum(f[k, j2, j1] for j2 in 0:instance.m if j2 == 0 || (j1 != j2 && instance.d_com[j1,j2] <= instance.rCom)))
            else # If connexity == :MultiFlow3hh

                # HAPS k1 and k2 can not communicate if they are not both deployed
                @constraint(model, [k1 in 1:instance.p, k2 in 1:instance.p; k1 != k2], q[k1, k2] <= sum(r[j, k1] for j in 1:instance.m))
                @constraint(model, [k1 in 1:instance.p, k2 in 1:instance.p; k1 != k2], q[k1, k2] <= sum(r[j, k2] for j in 1:instance.m))
                
                # If k1 is at position j1 and k2 is out of reach, q[k1, k2] = 0
                @constraint(model, [k1 in 1:instance.p, k2 in 1:instance.p, j1 in 1:instance.m; k1 != k2], q[k1, k2] <= 2 - r[j1, k1] - sum(r[j2, k2] for j2 in 1:instance.m if j2 != j1 && instance.d_com[j1, j2] > instance.rCom) - r[j1, k2])
                
                # No flow k if no HAPS in k
               @constraint(model, [k in 2:instance.p, k1 in 1:instance.p, k2 in 1:instance.p; k1 != k && k1 != k2], f[k, k1, k2] <= sum(r[j, k] for j in 1:instance.m))
                
                # The flow of HAPS k starts at HAPS 1 iff it is deployed
                @constraint(model, [k in 2:instance.p], sum(f[k, 1, k2] for k2 in 2:instance.p) == sum(r[j, k] for j in 1:instance.m))


                # y_j = sum_k r_j,k
                # relacher integrite des yj
                #sum_j r_jk <= 1
                # alpha_k <= alpha_k+1 ?
                # autre contrainte pour q_k1,k2 <= ... -> q_k1,k2 <= (1-r_j1k1)(1-r_j2k2) si j1 et j2 sont > rcom
                
                # No flow on (k1, k2) if the HAPS are not deployed or out of reach
                @constraint(model, [k in 2:instance.p, k1 in 1:instance.p, k2 in 1:instance.p; k1 != k && k1 != k2], f[k, k1, k2] <= q[k1, k2])

                # Flow conservation (source and think)
                # The flow of k reaching k is equal to that of flow k leaving 1
                @constraint(model, [k in 2:instance.p], sum(f[k, k1, k] for k1 in 1:instance.p if k1 != k) == sum(f[k, 1, k2] for k2 in 1:instance.p if k2 != 1))

                @constraint(model, [k in 2:instance.p, k1 in 1:instance.p; k1 != 1 && k1 != k], sum(f[k, k1, k2] for k2 in 1:instance.p if k2 != k1) == sum(f[k, k2, k1] for k2 in 1:instance.p if k2 != k1 && k2 != k ))
            end 
            
        else
            
            # If j is closed, it has no predecessor
            @constraint(model,  [j in 1:instance.m] , sum( x[j1,j] for j1 in 1:instance.m if j1 != j) <= y[j])
            
            # If j is closed, it has no successor
            # Remark: unlike the previous constraint there is no sum here as if y[j] = 1, j might have several successors
            @constraint(model,  [j in 1:instance.m, j1 in 1:instance.m ; j!=j1] , x[j,j1]<= y[j])
            
            # Number of connections = Number of opened sites - 1
            @constraint(model,  sum(x[j,j1] for j in 1:instance.m for j1 in 1:instance.m if xExist(j, j1, connexity)) == sum(y[j] for j in 1:instance.m) -1)
        end 

        if connexity == :MTZ

            # The order of the nodes is consistent with the arborescence
            @constraint(model,  [j in 1:instance.m, j1 in 1:instance.m ; j1!=j] , t[j1] >= t[j] + 1 - instance.p*(1-x[j,j1]))
        end 
    end


    # If a callback is used
    if connexity == :SubtoursCallback && !isRelaxation
        
        MOI.set(model, MOI.NumberOfThreads(), 1)

        cutsCount = 0
        lastNodeId = -1
    
        # Callback function
        function callbackSubTours(cb_data::CPLEX.CallbackContext, context_id::Clong)

            isInteger = isIntegerPoint(cb_data, context_id)
            isRelaxation = context_id == CPX_CALLBACKCONTEXT_RELAXATION

            # If the callback is called on an integer solution
            if isInteger || (isRelaxation && cutsCount <= maxCutsCount)

                valueP = Ref{Clong}()
                ret = CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODEUID, valueP)
                
                if lastNodeId != valueP[]
                    lastNodeId = valueP[]
                    cutsCount = 0
                end
                
                # Get the solution
                CPLEX.load_callback_variable_primal(cb_data, context_id)

                # Create the subtour problem
                yBar = callback_value.(cb_data, y)
                xBar = callback_value.(cb_data, x)
                
                sep, a, h = getSubPbModel(xBar, yBar, instance)
                set_silent(sep) 
                optimize!(sep)

                #Check if a violated constraint is found
                if (isInteger && objective_value(sep) > 0.001) || (isRelaxation && objective_value(sep) > 0.1)

                    # |links in S| <= |S| - 1
                    cstr = @build_constraint(sum(x[j,j1] for j in 1:instance.m-1 for j1 in j+1:instance.m if value(a[j]) > 1-1^-4 && value(a[j1])  > 1-1^-4)
                                             <= sum(y[j] for j in 1:instance.m if value(a[j]) > 1-1^-4) - sum(y[j] for j in 1:instance.m if value(h[j]) > 1 - 1^-4))

                    if isRelaxation
                        MOI.submit(model, MOI.UserCut(cb_data), cstr)
                        cutsCount += 1
                    else
                        MOI.submit(model, MOI.LazyConstraint(cb_data), cstr)
                    end  
                end
            end
        end

        # Add the callback to the model
        MOI.set(model, CPLEX.CallbackFunction(), callbackSubTours)

    end 

    return model, x, y, z, f, r
end


"""
Get the quadratic model of the localisation problem

Input
- instance: the instance solved
- isRelaxation: true if the variable are all continuous
- maxCutsCount : nombre maximal de coupes ajoutées à un noeud du branch-and-bound

Output:
- model: the model
- alpha, beta, z: its variables
"""
function getQuadraticModel(instance::Instance; time_limit::Int=-1, isRelaxation::Bool=false, threads::Int=-1)

    xMin = minimum(instance.clientPositions[:, 1])
    xMax = maximum(instance.clientPositions[:, 1])
    yMin = minimum(instance.clientPositions[:, 2])
    yMax = maximum(instance.clientPositions[:, 2])
    
    M = (xMax-xMin)^2+(yMax-yMin)^2

    ### Model declaration
    model = Model(CPLEX.Optimizer)
    #set_silent(model)
    if time_limit != -1 && time_limit > 5
        set_optimizer_attribute(model, "CPX_PARAM_TILIM", time_limit)
    end

    if threads != -1
        MOI.set(model, MOI.NumberOfThreads(), threads)
    end 
        
    ### Variables
    x = nothing
    z = nothing
    c = nothing
    
    if isRelaxation
        
        # 1 iff client i is covered by HAPS k
        @variable(model, 0 <= c[1:instance.n, 1:instance.p] <= 1)

        # 1 iff (j,j1) is a connection 
        @variable(model, 0 <= x[j in 1:instance.p, j1 in 1:instance.p ; j1!=j] <= 1)
        
    else         
        # 1 iff client i is covered by HAPS k
        @variable(model, c[1:instance.n, 1:instance.p], Bin)

        # 1 iff (j,j1) is a connection 
        @variable(model, x[j in 1:instance.p, j1 in 1:instance.p ; j1!=j], Bin)
        
    end
    
    # 1 iff client i is covered 
    @variable(model, 0 <= z[1:instance.n] <= 1)

    # Coordinates of the HAPS
    @variable(model, xMax >= alpha[k in 1:instance.p] >= xMin)
    @variable(model, yMax >= beta[k in 1:instance.p] >= yMin)
    
    # Index of site j in the connected arborescence
    @variable(model, t[j in 1:instance.p] >= 0)

    ### Objective
    @objective(model, Max, sum(z[i]  for i in 1:instance.n))

    ### Constraints

    # x_k,k2 = 0 if the distance between the HAPS is > rCom
    @constraint(model, [k in 1:instance.p, k2 in 1:instance.p; k != k2], (alpha[k] - alpha[k2])^2 + (beta[k]-beta[k2])^2 <= instance.rCom^2 * x[k, k2] + M * (1-x[k, k2]))
    
    # c_i,k = 0 if the distance between the HAPS and the client is > rCouv
    @constraint(model, [i in 1:instance.n, k in 1:instance.p], (alpha[k] - instance.clientPositions[i, 1])^2 + (beta[k]-instance.clientPositions[i, 2])^2 + instance.L^2 <= instance.rCouv^2 * c[i, k] + dFurthestHAPS(instance, i, xMin, xMax, yMin, yMax) * (1-c[i, k])) 
    
    # Client i is covered if at least one HAPS covers it
    @constraint(model, [i in 1:instance.n], z[i] <= sum(c[i, k] for k in 1:instance.p))

    # There is exactly p-1 variables x = 1
    @constraint(model, sum(x[k, k2] for k in 1:instance.p for k2 in 1:instance.p if k != k2) == instance.p-1)

    # A site does not have more than one input link
    @constraint(model,  [j in 1:instance.p] , sum(x[j1,j] for j1 in 1:instance.p if j1 != j) <= 1) 

    # The order of the nodes is consistent with the arborescence
    @constraint(model,  [j in 1:instance.p, j1 in 1:instance.p ; j1!=j] , t[j1] >= t[j] + 1 - instance.p*(1-x[j,j1])) 

    # Symmetry breaking
    @constraint(model, [k in 1:instance.p, k2 in k+1:instance.p], (yMax-yMin) * alpha[k] + beta[k] <= (yMax-yMin) * alpha[k2] + beta[k2])

    return model, alpha, beta, z, c, t, x
end 

"""
Create a linear relaxation of the model.
Generate subtour inequalities.
Set the variables to binary and then set add the MTZ variables and constraints.

Input:
- instance: the instance to solve.
- (optional) time_limit: resolution time limit (-1 if no limit) (default -1)
- (optional) maxCuts: maximal number of subtour cuts added at the root (used only if connexity = is :RelaxSubtoursThenMTZ) (default 50)

Output:
- model: the model
- x, y, z: its variables
"""
function generateRelaxedSubtours(instance::Instance, maxCuts::Int, time_limit::Int; isRelaxation::Bool=false, upperBound::Int64=typemax(Int64), threads::Int=-1, ppResults=nothing)
    
    model, x, y, z = getModel(instance, connexity = :SubtoursCallback, isRelaxation = true, ppResults = ppResults)

    if threads != -1
        MOI.set(model, MOI.NumberOfThreads(), threads)
    end 
    
    results = Dict{String, Any}()
    
    sSets = Array{Array{Int}}([])
    best_lower_bound=0
    isOptimal =  abs(best_lower_bound - upperBound) < 1E-4
    nbCuts = 0
    startingTime = time()
    remainingTime = time_limit

    while !isOptimal && nbCuts < maxCuts && (remainingTime > 5 || time_limit == -1)

        if time_limit != -1 
            set_optimizer_attribute(model, "CPX_PARAM_TILIM", remainingTime)
        end
        
        set_silent(model)
        optimize!(model)     
        isOptimal = true
        
        if termination_status(model) == MOI.OPTIMAL
            best_lower_bound=objective_value(model)

            #Separate subtour elimination constraints
            sep, a, h = getSubPbModel(value.(x), value.(y), instance)
            set_silent(sep) 
            optimize!(sep)
            
            exist_violated_cut= false

            #Check if some subtour constraint is not verified
            if objective_value(sep) > 0.001
                exist_violated_cut = true

                S = Vector{Int}()

                j0 = -1
                for j in 1:instance.m
                    if value(a[j]) > 1 - 1^-4
                        push!(S, j)
                    end
                    if value(h[j]) > 1 - 1^-4
                        j0 = j
                    end
                end

                push!(S, j0)
                push!(sSets, S)
                
                # |links in S| <= |S| - 1
                @constraint(model, sum(x[S[j],S[j1]] for j in 1:length(S)-1 for j1 in j+1:length(S)-1)
                            <= sum(y[S[j]] for j in 1:length(S)-1) - y[S[end]])
                nbCuts += 1
            end
            isOptimal = !exist_violated_cut || abs(best_lower_bound - upperBound) < 1E-4
        end

        elapsedTime = time() - startingTime
        remainingTime = time_limit - elapsedTime
    end
    results["rootRelaxationAfterST"] = best_lower_bound
    results["stCutsRootCount"] = nbCuts
    results["rootSTSeparationTime"] = time()-startingTime
    results["isRootOptimal"] = isOptimal

    if !isRelaxation && (remainingTime > 5 || time_limit == -1)

        intTimeLimit = -1

        if time_limit != -1
            intTimeLimit = round(Int, remainingTime)
        end
        
        model, x, y, z = getModel(instance, connexity = :MTZ, time_limit = remainingTime, ppResults = ppResults) 

        # Add all the generated inequalities
        for S in sSets 
            @constraint(model, sum(x[S[j],S[j1]] + x[S[j1],S[j]] for j in 1:length(S)-1 for j1 in j+1:length(S)-1)
                        <= sum(y[S[j]] for j in 1:length(S)-1) - y[S[end]])
        end
    end 
    #unset_silent(model)

    return model, x, y, z, results
    
end 

"""
Return the subproblem model associated with an instance and values of variables x and y.

Input:
- vx, vy: value of x and y
- instance: the instance to solve

Output:
- sep: the model
- a: its variables
"""
function getSubPbModel(vx, vy, instance)

    sep = Model(CPLEX.Optimizer) 

    # j in S?
    @variable(sep, a[1:instance.m], Bin)

    # j and j1 in S? (= a[j] * a[j1])
    @variable(sep, w[j in 1:instance.m, j1 in 1:instance.m; j<j1] >= 0)

    # Linearization of w_j,j1 = a_j * a_j1
    @constraint(sep, inf1[j=1:instance.m-1,j1=j+1:instance.m], w[j,j1] <= a[j])
    @constraint(sep, inf2[j=1:instance.m-1,j1=j+1:instance.m], w[j,j1] <= a[j1])

    # h is the selected open site in S?
    @variable(sep, h[1:instance.m], Bin)

    # Maximize the number of links in S
    # - # open sites in S
    # + exists open sites in S?
    @objective(sep, Max,  sum(value(vx[j,j1])*w[j,j1] for j in 1:instance.m-1 for j1 in j+1:instance.m) - sum(value(vy[j])*a[j] for j in 1:instance.m) + sum(value(vy[j])*h[j] for j in 1:instance.m)  )
    @constraint(sep, inf3[j=1:instance.m], h[j] <= a[j])
    @constraint(sep, sum(h[j] for j in 1:instance.m) ==1)
    
    return sep, a, h
end 

"""
Determines if a callback is called because an integer solution is found

Intput:
- cb_data, context_id: callback context

Output:
- true if the callback is called because CPLEX has found an integer point
"""
function isIntegerPoint(cb_data::CPLEX.CallbackContext, context_id::Clong)

    # context_id  == CPX_CALLBACKCONTEXT_CANDIDATE si le  callback est
    # appelé dans un des deux cas suivants :
    # cas 1 - une solution entière a été obtenue; ou
    # cas 2 - une relaxation non bornée a été obtenue
    if context_id != CPX_CALLBACKCONTEXT_CANDIDATE
        return false
    end

    # Pour déterminer si on est dans le cas 1 ou 2, on essaie de récupérer la
    # solution entière courante
    ispoint_p = Ref{Cint}()
    ret = CPXcallbackcandidateispoint(cb_data, ispoint_p)

    # S'il n'y a pas de solution entière
    if ret != 0 || ispoint_p[] == 0
        return false
    else
        return true
    end
end

function solveGreedy(instancePath::String; isConnected::Bool=true, isStochastic::Bool=true, p::Int = -1, time_limit::Int=-1)

    instance = Instance(instancePath)
    return solveGreedy(instance, isConnected = isConnected, isStochastic = isStochastic, p = p, time_limit = time_limit)
end 

"""
Greedily open centers which cover the most clients

Input
- instance: the instance to solve
- (optional) isConnected: true if the drones must be able to communicate with each others; false otherwise
"""
function solveGreedy(instance::Instance; isConnected::Bool=true, isStochastic::Bool=true, p::Int = -1, time_limit::Int=-1)

    if p != -1
        instance.p = p
    end
    
    startingTime = time()

    # isClientCovered[i] = 1 if client i is covered in the greedy solution; 0 otherwise
    isClientCovered = Vector{Float64}(zeros(instance.n))
    
    # isSiteOpened[j] = 1 if site j is opened in the greedy solution; 0 otherwise
    isSiteOpened = Vector{Float64}(zeros(instance.m))

    # isSiteReachble[i] is true iff site i is not currently in the solution and can reach one of the site in the current solution
    # (initially there is no site in the solution)
    isSiteReachable = Vector{Bool}(zeros(instance.m))
    
    # isClientCoverable[i, j] = 1 if client i is not already covered by a site and is reachable from drone j; 0 otherwise
    isClientCoverable = Array{Int, 2}(instance.d .<= instance.rCouv)

    # coverableClientCount[i] = number of not currently covered clients that can be covered by site i
    coverableClientsCount =  vec(sum(isClientCoverable, dims = 1))
    
    # Get a site which covers the maximal number of clients
    connectedClientsCount = 0
    bestSite = -1

    if isStochastic
        # Select a site randomly
        # (the more new clients a site enables to cover, the more probable it is to be selected)
        bestSite = sample(Weights(exp.(coverableClientsCount))) 
        coverableClientsCount = coverableClientsCount[bestSite]
    else
        # Select a site which enables to cover the most clients
        connectedClientsCount, bestSite = findmax(coverableClientsCount)
    end

    isSiteOpened[bestSite] = 1

    # Add it to the opened sites
    openedSites = Vector{Int}([bestSite])     

    # Remove the covered clients 
    for i in 1:instance.n
        if isClientCoverable[i, bestSite] == 1
            isClientCoverable[i, :] .= 0
            isClientCovered[i] = 1
        end
    end

    # Update the coverable client count of each site
    coverableClientsCount =  vec(sum(isClientCoverable, dims = 1))
    
    # Update the list of sites that can be reached
    for j in 1:instance.m
        if j != bestSite && (!isConnected || instance.d_com[bestSite, j] <= instance.rCom)
            isSiteReachable[j] = true
        end 
    end
    
    emptyVector = Vector{Int}(zeros(instance.m))

    elapsedTime = time() - startingTime
    remainingTime = time_limit - elapsedTime
    
    # While sites can be added to the solution
    while length(openedSites) < instance.p && isSiteReachable != emptyVector && (remainingTime > 0 || time_limit == -1)

        # Select a site  
        if isStochastic
            # Randomly get a site which is reachable and not already imposed to be opened
            # (the more new clients it enables to cover, the more likely it is to be picked)
            # (0.1 is added to avoid having a weight equal to 0 for all sites which could lead to selecting an unreachable site)
            bestSite = sample(Weights((0.1 .+ exp.(coverableClientsCount)) .* isSiteReachable))
            bestClientCount = coverableClientsCount[bestSite] 
        else
            # Find the best reachable site (i.e., the site that covers the most of uncovered clients) 
            # (0.1 is added to avoid having a weight equal to 0 for all sites which could lead to selecting an unreachable site)
            bestClientCount, bestSite = findmax((coverableClientsCount .+ 0.1) .* isSiteReachable .* (1 .- isSiteOpened))
        end

        connectedClientsCount += bestClientCount
        # Add it
        push!(openedSites, bestSite)
        isSiteOpened[bestSite] = 1
        isSiteReachable[bestSite] = false

        # Remove the covered clients 
        for i in 1:instance.n
            if isClientCoverable[i, bestSite] == 1
                isClientCoverable[i, :] .= 0
                isClientCovered[i] = 1
            end
        end

        # Update the coverable client count of each site
        coverableClientsCount =  vec(sum(isClientCoverable, dims = 1))
        
        # Update the list of sites that can be reached
        for j in 1:instance.m
            if isSiteOpened[j] == 0 && (!isConnected || instance.d_com[bestSite, j] <= instance.rCom)
                isSiteReachable[j] = true
            end 
        end
        
        elapsedTime = time() - startingTime
        remainingTime = time_limit - remainingTime
        
    end
    
    results = Dict{String, Any}()
    results["resolutionTime"] = time() - startingTime
    results["isSiteOpened"] = isSiteOpened
    results["isClientCovered"] = isClientCovered
    results["objective"] = sum(isClientCovered)
    results["p"] = instance.p

    return results
    
end 


"""
Greedily open centers which cover the most clients.
The site opened at a given iteration is:
- reachable by at least one already opened site (except for the first iteration);
- chosen randomly with weights equal to the value of the linear relaxation of variables y.

Input
- instancePath: path of the instance
- (optional) isConnected: true if the drones must be able to communicate with each others; false otherwise
"""
function solveGreedyRelax(instancePath::String; isConnected::Bool=true, isStochastic::Bool=true, p::Int=-1, time_limit::Int=-1)

    instance = Instance(instancePath)
    
    if p != -1
        instance.p = p
    end 

    connexity = :MTZ

    if !isConnected
        connexity = :None
    end

    if time_limit != -1 && time_limit < 5
        time_limit = -1
    end
    
    startingTime = time()
    
    # isSiteOpened[j] = 1 if site j is opened in the greedy solution; 0 otherwise
    isSiteOpened = Vector{Float64}(zeros(instance.m))


    # isSiteReachble[i] is true iff site i is not currently in the solution and can reach one of the site in the current solution
    # (initially there is no site in the solution)
    isSiteReachable = Vector{Bool}(zeros(instance.m))

    # canSiteBeOpened[j] is true unless site j is forced to be open or it has previously been opened and it had lead to an infeasibility
    canSiteBeOpened = Vector{Bool}(ones(instance.m))

    m, x, y, z = getModel(instance, connexity = connexity, isRelaxation = true; time_limit = time_limit)
    set_silent(m)
    optimize!(m)
    
    # Get a site which covers the maximal number of clients
    bestSite = -1

    if isStochastic
        # Randomly get a site (only if the site is not already imposed to be opened or lead to an infeasibility)
        bestSite = sample(Weights(exp.(1 .+ 10 .* JuMP.value.(y)))) 
    else
        # Find the best site (i.e., the site for which y[j] is maximal) 
        bestValue, bestSite = findmax(JuMP.value.(y))
    end

    isSiteOpened[bestSite] = 1
    canSiteBeOpened[bestSite] = false
    @constraint(m, y[bestSite] == 1)
        
    # At each iteration, a site is opened
    iterationCount = 1
    openedSites = 1

    # Update the list of sites that can be reached
    for j in 1:instance.m
        if j != bestSite && (!isConnected || instance.d_com[bestSite, j] <= instance.rCom)
            isSiteReachable[j] = true
        end 
    end

    emptyVector = Vector{Int}(zeros(instance.m))
    candidateSites = canSiteBeOpened .* isSiteReachable

    elapsedTime = time() - startingTime
    remainingTime = time_limit - elapsedTime
    
    # While sites can be opened (<p)
    # and there exists connected sites to open
    # and the time is not limited or there is at least 5 seconds remaining
    while openedSites < instance.p && candidateSites != emptyVector && (time_limit == -1 || remainingTime >= 5)

        # Reduce the remaining time 
        if time_limit != -1
            set_optimizer_attribute(m, "CPX_PARAM_TILIM", remainingTime)
        end
        
        optimize!(m)

        if primal_status(m) == MOI.FEASIBLE_POINT 

            if isStochastic
                # Randomly get a site (only if the site is not already imposed to be opened or lead to an infeasibility)
                # (0.1 is added to avoid opening a site which can not be opened)
                bestSite = sample(Weights((0.001 .+ exp.(1 .+ 10 .* (JuMP.value.(y)))) .* candidateSites)) 
            else
                # Find one of the best site (i.e., one of the sites for which y[j] is maximal) 
                bestValue, bestSite = findmax((0.001 .+ JuMP.value.(y)) .* candidateSites)
            end
            
            openedSites += 1
            isSiteOpened[bestSite] = 1
            canSiteBeOpened[bestSite] = false
            isSiteReachable[bestSite] = false

            @constraint(m, y[bestSite] == 1)
            
            # Update the list of sites that can be reached
            for j in 1:instance.m
                if isSiteOpened[j] == 0 && (!isConnected || instance.d_com[bestSite, j] <= instance.rCom)
                    isSiteReachable[j] = true
                end 
            end

            iterationCount += 1
            candidateSites = canSiteBeOpened .* isSiteReachable
        end
        
        elapsedTime = time() - startingTime 
        remainingTime = time_limit - elapsedTime  
    end

    # isClientCovered[i] = 1 if client i is covered in the greedy solution; 0 otherwise
    isClientCovered = Vector{Float64}(zeros(instance.n))

    for i in 1:instance.n
        isCovered = false

        j = 1
        while j < instance.m && !isCovered
            if isSiteOpened[j] == 1 && instance.d[i,j] <= instance.rCouv
                isCovered = true
            end
            j += 1
        end

        isClientCovered[i] = isCovered
    end 

    results = Dict{String, Any}()
    results["resolutionTime"] = time() - startingTime
    results["isSiteOpened"] = isSiteOpened
    results["isClientCovered"] = isClientCovered
    results["objective"] = sum(isClientCovered)
    results["p"] = instance.p

    return results
end 

"""
True if a variable x_j1,j2 exists in a model
"""
function xExist(j1::Int, j2::Int, connexity::Symbol)
    if connexity == :SubtoursCallback
        return j1 < j2
    else
        return j1 != j2
    end
end 

function doPreprocessing(instance)

    # isClientCoverable[i, j] = 1 if client i is coverable by site j
    isClientCoverable = Array{Bool, 2}(instance.d .<= instance.rCouv)

    # isCommunicationPossible[i, j] = true iff site i and j can communicate
    isCommunicationPossible = Array{Bool, 2}(instance.d_com .<= instance.rCom)

    # coverableClientsCount[j] = number of clients covered by position j
    coverableClientsCount =  vec(sum(isClientCoverable, dims = 1))


    ## Find for each site, a site which covers the same clients

    isDominated = Vector{Bool}(zeros(instance.m))
    

    # A line [j1 j2] of priorities means that site j2 can only be used if site j1 is also used
    priorities = Matrix{Int}(zeros(0, 2))

    # A line [j1 j2 s] of communicationPriorities means that site j2 can only be used it communicates with a site unreachable by j1 (all these sites are contained in the Vector{Int} s)
    communicationPriorities = Vector{Dict{String, Any}}()
    
    for j1 in 1:instance.m
        for j2 in j1+1:instance.m

            if !isDominated[j1] && !isDominated[j2]
                diffClients = isClientCoverable[:, j1] .- isClientCoverable[:, j2]
                j1CoversJ2Clients = minimum(diffClients) >= 0
                j2CoversJ1Clients = maximum(diffClients) <= 0

                diffSites = isCommunicationPossible[j1, :] .- isCommunicationPossible[j2, :]
                j1CoversJ2Sites = minimum(diffSites) >= 0
                j2CoversJ1Sites = maximum(diffSites) <= 0
                
                if j1CoversJ2Clients
                    if j1CoversJ2Sites # If j1 covers all clients and sites covered by j2
                        isDominated[j2] = true
                        println(j1, " >= ", j2)
                        
                    elseif j2CoversJ1Clients && j2CoversJ1Sites # If j2 covers the same clients than j1 and all its sites
                        isDominated[j1] = true
                        println(j2, " >= ", j1)
                    else

                        # If j1:
                        # - covers all clients covered by j2; and
                        # - does not cover all sites of j2; and
                        # - covers strictly more clients than j2 or has not all its sites covered by j2
                        # Then j2 should only be used if it communicates with a site unreachable by j1
                        unreachableSites = vec(findall(isCommunicationPossible[j1, :] - isCommunicationPossible[j2, :] .== -1))
                        d = Dict{String, Any}()
                        d["site1"] = j1
                        d["site2"] = j2
                        d["unreachableSites"] = unreachableSites
                        push!(communicationPriorities, d)
                    end 
                    
                elseif j2CoversJ1Clients
                    if j2CoversJ1Sites # If j2 covers all clients and sites covered by j1
                        isDominated[j1] = true
                        println(j2, " >= ", j1)
                    else
                        # If j2:
                        # - covers more clients than j1; and
                        # - does not cover all sites of j1
                        # Then j1 should only be used if it communicates with a site unreachable by j1
                        unreachableSites = vec(findall(isCommunicationPossible[j2, :] - isCommunicationPossible[j1, :] .== -1))
                        d = Dict{String, Any}()
                        d["site1"] = j2
                        d["site2"] = j1
                        d["unreachableSites"] = unreachableSites
                        push!(communicationPriorities, d)
                    end
                else # If the clients covered by j1 are not all covered by j2 an vice versa

                    ## !! Not always true
                    #= # If j1 covers all sites covered by j2 and more clients
                    if j1CoversJ2Sites && sum(isClientCoverable[:, j1]) >= sum(isClientCoverable[:, j2])
                        # Then j2 can only by used if j1 is also used
                        priorities = vcat(priorities, [j1 j2])
                    elseif j2CoversJ1Sites && sum(isClientCoverable[:, j1]) <= sum(isClientCoverable[:, j2])
                        priorities = vcat(priorities, [j2 j1])
                    end
                    =#
                end
            end
        end
    end

    results = Dict{String, Any}()
    results["priorities"] = priorities
    results["comPriorities"] = communicationPriorities
    results["dominated"] = findall(isDominated .== true)
    return results
end 

"""
Returns the furthest distance betwee client i and an HAPS
(i.e., with client i and a corner of the map)
"""
function dFurthestHAPS(instance::Instance, i::Int, xMin, xMax, yMin, yMax)
    x = instance.clientPositions[i, 1]
    y = instance.clientPositions[i, 2]
    
    return maximum([(x-xMin)^2+(y-yMin)^2, (x-xMax)^2+(y-yMin)^2, (x-xMin)^2+(y-yMax)^2, (x-xMax)^2+(y-yMax)^2])
    
end 

function isFVariable(s, t, j1, j2, instance)
    return instance.d_com[j1,j2] <= instance.rCom && j1 != j2 && s < t && j1 != t && j2 != s
end

"""
Iteratively solve the relaxation of a formulation and separate clique inequalities to tighten its value

Input:
- redundancyRatio: percentage in [0, 1] of units in a new clique that must not have appeared in any previous clique
  = 0 : no constraint
  = 0.5 : half the units in a clique must not have appeared in a previous clique
  = 1 : all the units must not have appeared in a previous clique
"""
function tightenRelaxationByCliques(instance::Instance; connexity::Symbol=:MTZ, p::Int=-1, redundancyRatio::Float64=Float64(0.0), time_limit::Int=-1)

    startingTime = time()
    if p != -1
        instance.p = p
    end

    # isInAClique[u] is true if u is in at least one of the clique found
    # (only used if the number of new units in a clique is constrained)
    isInAClique = Vector{Bool}([])

    if redundancyRatio > 0
        isInAClique = Vector{Bool}(zeros(instance.n))
    end 
    
    model = nothing
    z = nothing

    if connexity == :QuadMTZ
        model, alpha, beta, z, c, t, x =  getQuadraticModel(instance, isRelaxation = true)
    else 
        model, x, y, z, f, r = getModel(instance, connexity = connexity, isRelaxation = true)
    end 
    set_silent(model) 

    violatedCliqueFound = true

    cliques = Vector{Vector{Int}}()

    areClientsClose = Matrix{Bool}(zeros(instance.n, instance.n))

    closeThreshold = 2 * sqrt(instance.rCouv^2 - L^2)
    for u1 in 1:instance.n
        for u2 in u1+1:instance.n
            dist = sqrt((instance.clientPositions[u1, 1] - instance.clientPositions[u2, 1])^2+(instance.clientPositions[u1, 2] - instance.clientPositions[u2, 2])^2)
            if dist <= closeThreshold
                areClientsClose[u1, u2] = true
                areClientsClose[u2, u1] = true
            end 
        end
    end

    firstRelaxation = -1
    lastRelaxation = -1

    # Number of iteration at which the improvement of the relaxation is tested
    testFrequency = 100

    # True if the relaxation has been improved sufficiently at the last test
    sufficientImprovement = true
    
    # Value of the relaxation when the improvement was last tested
    lastTestedRelaxation = instance.n

    iterationCount = 0

    elapsedTime = time() - startingTime

    while violatedCliqueFound && sufficientImprovement && (time_limit == -1 || elapsedTime <= time_limit/2)
        iterationCount += 1
        optimize!(model)

        lastRelaxation = JuMP.objective_value(model)
        if firstRelaxation == -1
            firstRelaxation = lastRelaxation
        end 
        println("Bound: ", round(lastRelaxation, digits = 2))
        
        newClique = cliqueSeparation(instance, JuMP.value.(z), areClientsClose, isUnitForbidden=isInAClique, redundancyRatio=redundancyRatio)

        if length(newClique) > 0
            println("Found clique : ", newClique)
            push!(cliques, newClique)
            @constraint(model, sum(z[u] for u in newClique) <= instance.p)

            if redundancyRatio > 0
	       for unit in newClique
                  isInAClique[unit] = true
	       end
            end 
        else
            println("No clique found")
            violatedCliqueFound = false
        end

        if rem(iterationCount, testFrequency) == 0
            if (lastTestedRelaxation - lastRelaxation)/lastRelaxation <= 0.01
                sufficientImprovement = false
            end
            lastTestedRelaxation = lastRelaxation
        end

        elapsedTime = time() - startingTime
    end

    return cliques, firstRelaxation, lastRelaxation
end 

"""
Let G=(V, E, z) be the graph such that:
- V are the clients
- (i,j) is in E if d_i,j > rCouv
- z_i is the weight of client i (1 iff i is covered)

Find a clique of G of maximal weight

Input
- instance: the instance
- z: the z values
- areClientsClose[u1, u2] is true if the distance between clients u1 and u2 is <= 2*sqrt(rCouv^2-L^2)
- (optional) isUnitForbidden: vector of size instance.n such that isUnitForbidden[u] is true if unit u is ignored in the search of clique (e.g., to ensure that the units in the clique have not already been included in a previous clique).
- (optional) redundancyRatio: percentage in [0, 1] of units in a new clique that must not have appeared in any previous clique
  = 0 : no constraint
  = 0.5 : half the units in a clique must not have appeared in a previous clique
  = 1 : all the units must not have appeared in a previous clique

Output
- a clique of size >p if any exists; an empty vector otherwise
"""
function cliqueSeparation(instance::Instance, z::Vector{Float64}, areClientsClose::Matrix{Bool}; isUnitForbidden::Vector{Bool}=Vector{Bool}([]), redundancyRatio::Float64=Float64(0.0))

    cm = Model(CPLEX.Optimizer)
    set_silent(cm) 

    # nu[i] = 1 iff client i is in the clique
    @variable(cm, nu[1:instance.n ], Bin)

    @objective(cm, Max, sum(z[u] * nu[u] for u in 1:instance.n))
    @constraint(cm, [u1 in 1:instance.n, u2 in u1+1:instance.n; areClientsClose[u1, u2]], nu[u1]+ nu[u2] <= 1)

    if length(isUnitForbidden) > 0

        # The number of new units in the clique must be greater than redundancyRatio% of the units in the clique
        # (here "new" = have not already been included in a previously found clique)
        @constraint(cm, sum(nu[u] for u in 1:instance.n if !isUnitForbidden[u]) >= redundancyRatio * sum(nu[u] for u in 1:instance.n))
    end
    
    optimize!(cm)

    clique = Vector{Int}([])

    if primal_status(cm) == MOI.FEASIBLE_POINT && JuMP.objective_value(cm) > instance.p
        for u in 1:instance.n
            if JuMP.value(nu[u]) > 0.9
                push!(clique, u)
            end
        end 
    end

    return clique

end 
