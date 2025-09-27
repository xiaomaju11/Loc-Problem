"""
Add to a model a clique inequality that corresponds to a set of units
"""
function getCliqueInequality(model::Model, instance::Instance, clique::Vector{Int64}, dVariables::Dict{Symbol, Any})
    return getCliqueInequality(InequalityFamily(), model, instance, clique, dVariables)
end
    

"""
Add to a model a clique inequality that corresponds to a set of units
"""
function getCliqueInequality(iFamily::InequalityFamily, model::Model, instance::Instance, clique::Vector{Int64}, dVariables::Dict{Symbol, Any})
    if haskey(dVariables, :z)
        return @build_constraint(sum(dVariables[:z][u] for u in clique) <= instance.p)
    else
        return @build_constraint(sum(dVariables[:z2][j, u] for u in clique, j in 1:instance.m) <= instance.p)
    end
end 

"""
Try to find an incompatible clusters inequality that separate a fractional solution
- dVariables: dictionary containing the variables used. It must contain one of the two next entries:
  - :z for single indexed z variables (z[u]=1 iff u is covered); or 
  - :z2 for doubled indexed z variables (z[j, u]=1 iff u is covered by site j); and it can additionally contain 
  - (optional) :y for the sites variables (y[j]=1 iff site j has an HAPS)

"""
function cliqueSeparation(iFamily::InequalityFamily, instance::Instance, model::Model, dVariables::Dict{Symbol, Any}; cb_data=nothing)
    
    if size(instance.areClientsClose, 1) == 0
        computeCloseClients!(instance)
    end

    #@show instance.areClientsClose
    
    # Set the weight used to pick the order in which the units are considered
    # (the greater the value of z, the more likely the unit is to be picked)
    variablesValues = nothing

    if haskey(dVariables, :z)
        variablesValues = getValue.(dVariables[:z], cb_data)
    else
        variablesValues = Vector{Float64}(ones(instance.n))
        for u in 1:instance.n
            variablesValues[u] = getValue.(dVariables[:z][:, u], cb_data)
        end 
    end

    #@show round.(variablesValues, digits = 2)
        
    # isInAClique[u] is true if client u is in one of the found cliques
    isInAClique = Vector{Bool}(zeros(instance.n))

    # List of all the cliques found
    cliques = Vector{Vector{Int}}([])

    # Violation of each clique found
    cliquesViolation = Vector{Float64}([])

    # Only units with positive value of the variables are considered as seeds
    # seed = first element added to a clique
    seedWeights = copy(variablesValues)

    # Increase all weights to allow units with variables equal to 0 in the cliques (even if they cannot be seeds)
    initialCliqueWeights = copy(variablesValues) .+ 1E-3
    
    remainingSeedCandidates = count(x -> x>1E-4, variablesValues)
    
    maxCutsAdded = min(iFamily.maxCutsAddedInOneIteration, iFamily.maxCutsAddedInNextIteration)

    maxCutsFound = 3 * maxCutsAdded

    time_limit = -1

    if haskey(iFamily.additionalData, :time_limit)
        time_limit = iFamily.additionalData[:time_limit]
    end

    maxSeparationTimeWithCut = -1
    maxSeparationTimeWithoutCut = -1

    if haskey(iFamily.additionalData, :relaxationTime)        
        maxSeparationTimeWithCut = 10 * iFamily.additionalData[:relaxationTime]
        maxSeparationTimeWithoutCut = 50 * iFamily.additionalData[:relaxationTime]
    end 
    
    startingTime = time()
    elapsedTime = 0

    #@show remainingSeedCandidates
    
    # While:
    # - the maximal number of cut found is not reached; and
    # - all the elements with non null variable value have not been used as seeds; and
    # - the separation time is lower than the relaxation time or no cut is found; and
    # - the separation time is lower than 5*(the relaxation time); and
    # - the time limit is not reached.
    while length(cliques) < maxCutsFound && remainingSeedCandidates > 0 && (maxSeparationTimeWithCut == -1 ||  elapsedTime < maxSeparationTimeWithCut || length(cutsFound) == 0) && (maxSeparationTimeWithoutCut == -1 || elapsedTime < maxSeparationTimeWithoutCut) && (time_limit == -1 || elapsedTime <= time_limit)

        # Start a clique with a random unit
        seedUnit = sample(Weights(seedWeights)) 
        clique = Vector{Int}([seedUnit])

        currentCliqueWeights = copy(initialCliqueWeights)
        remainingCliqueCandidates = instance.n
            
        for u in 1:instance.n
            if instance.areClientsClose[u, seedUnit] || u == seedUnit#### belong to the same clique
                currentCliqueWeights[u] = 0
                remainingCliqueCandidates -= 1
                #@show u
            end
        end 

        #@show seedUnit, remainingCliqueCandidates
        #@show round.(currentCliqueWeights, digits = 2)

        while remainingCliqueCandidates > 0

            # Select a unit to add to the clique
            nextUnit = sample(Weights(currentCliqueWeights))
            
            # Since only valid candidates have a weight >0, nextUnit can be added to the clique
            push!(clique, nextUnit)

            if currentCliqueWeights[nextUnit] == 0
                println("Error: unit of weight 0 added to the clique (unit id: ", nextUnit, ")")
                readline()
            end

            for u in 1:instance.n
                if (instance.areClientsClose[u, nextUnit] || u == nextUnit) && currentCliqueWeights[u] > 0
                    currentCliqueWeights[u] = 0
                    remainingCliqueCandidates -= 1
                end 
            end 
            #@show nextUnit, remainingCliqueCandidates
        #@show round.(currentCliqueWeights, digits = 2)
        end

        cutLHS = 0.0
        cliqueId = 1

        while cutLHS <= instance.p + 1E-4 && cliqueId <= length(clique)
            cutLHS += variablesValues[clique[cliqueId]]###not cliques?
            cliqueId += 1
        end

        if cutLHS > instance.p + 1E-4
            push!(cliques, clique)
            push!(cliquesViolation, cutLHS)

            for unit in clique
                seedWeights[unit] = 0
                remainingSeedCandidates -= 1
            end

            #println("Adding clique: ", clique)
        else
            seedWeights[seedUnit] = 0
            remainingSeedCandidates -= 1 
        end
        
        elapsedTime = time() - startingTime 
    end # while length(cliques) < maxCutsAdded && remainingSeedCandidates > 0

    cliques = cliques[sortperm(cliquesViolation, rev = true)]

    #println("+" , length(cliques), " clique inequalities: ", length.(cliques))
    #@show cliques
    #readline()
    return cliques[1:min(end, maxCutsAdded)]
end 
