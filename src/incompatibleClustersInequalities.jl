"""
Add incompatible clusters inequalities to the formulation.
Let {C_i} be a set of clusters of clients (it can also contain site positions) such that two elements in different clusters cannot both be in a solution (they are too far away). The corresponding incompatible cluster inequality create one variable zC_i for each cluster C_i the constraints: 
sum_i zC_i <= 1
z_u <= zC_i for all unit in C_i
y_j <= zC_i for all site in C_i

Input:
- instance: the instance
- model: the model
- dVariables: dictionary containing the variables used. It must contain one of the two following entries:
  - :z for single indexed z variables (z[u]=1 iff u is covered); or 
  - :z2 for doubled indexed z variables (z[j, u]=1 iff u is covered by site j); and it can additionally contain 
  - :y for the sites variables (y[j]=1 iff site j has an HAPS)
"""
function generateAllCouplesICConstraints(instance::Instance, model::Model, dVariables::Dict{Symbol, Any})

    results = Dict{String, Any}()
    results["ciIneqCount"] = 0
    startingTime = time()

    y, z, z2 = getVariables(dVariables)
    positionCount = instance.n

    if y != nothing
        positionCount = instance.n + instance.m
    end 
    
    # Ajouter une variable zC par cluster
    # zC >= zc pour tout c du cluster
    # Dans la contrainte on a sum_C zC <= 1

    # isInAConstraint[u] is true if unit u is in one of the constraints previously found
    isInAConstraint = Vector{Bool}(zeros(positionCount))

    areInDifferentClusterOfAConstraint = Matrix{Bool}(zeros(positionCount, positionCount))

    # Dictionary that will be filled as new variables representing clusters are added to the problem
    # (used to avoid creating several variables representing the exact same cluster)
    previousClustersDictionary = Dict{Symbol, Any}()

    areIncompatible = getIncompatibilities(instance, positionCount)

    if areIncompatible != nothing
        
        # For each unit or site
        for u1 in 1:positionCount
            for u2 in u1+1:positionCount

                if !areInDifferentClusterOfAConstraint[u1, u2] && areIncompatible[u1, u2]

                    # Create the first two sets with u1 and u2
                    initialSets = Vector{Vector{Int}}()
                    push!(initialSets, Vector{Int}([u1]))
                    push!(initialSets, Vector{Int}([u2]))                    

                    sets = generateIC(instance, areIncompatible, initialSets=initialSets)

                    # Set which couple of elements have been in different sets
                    for setId1 in 1:length(sets)
                        set1 = sets[setId1]
                        for setId2 in setId1+1:length(sets)
                            set2 = sets[setId2]

                            for element1 in set1
                                for element2 in set2
                                    areInDifferentClusterOfAConstraint[element1][element2] = true
                                    areInDifferentClusterOfAConstraint[element2][element1] = true
                                end
                            end 
                        end
                    end

                    # Add the corresponding inequality and update the previous clusters dictionary if new cluster variables are added
                    addICInequality(model, sets, dVariables, additionalData=previousClustersDictionary)
                    results["ciIneqCount"] += 1
                    
                end # if !areInDifferentClusterOfAConstraint[u1, u2]
            end # for u2 in u1+1:positionCount
        end # for u1 in 1:positionCount
    end # if maximalDistance < sqrt((xMax-xMin)^2+(yMax-yMin)^2)

    results["ciTime"] = time() - startingTime
    results["ciNewVarCount"] = length(previousClustersVariables)
    return results
end 



"""
Generate an incompatible cluster constraint by randomly selecting an element to add in a new set or in an existing set.
The elements are randomly picked with a weight equal to initialElementsWeight.

Input
- instance: the instance
- areIncompatible[i, j] is true if elements i and j (either units and/or HAPS sites) cannot be in the same solution. If an index is <= n it corresponds to a unit, otherwise, it corresponds to a site.
- (optional) initialElementsWeights: weight of each element (the greater the weight, the more likely the element will be considered to be added to the inequality). Uniform weight are consiered if this argument is not specified
- initialSets: initial sets if some elements are initially imposed to be in a set (enables to ensure that one or several elements are in the constraint)
"""
function generateIC(instance::Instance, areIncompatible::Matrix{Bool}; initialElementsWeight::Vector{Float64}=Vector{Int}([]), initialSets::Vector{Vector{Int}}=Vector{Vector{Int}}([]))

    positionCount = size(areIncompatible, 1)

    # If no initial element weights are provided
    if length(initialElementsWeight) == 0

        # All elements have the same weight
        initialElementsWeight = ones(positionCount)
    end
    
    sets = initialSets

    # hasCompatibilityWithSet[s][u] is true if unit u is compatible with at least one unit in set s
    hasCompatibilityWithSet = Vector{Vector{Bool}}()

    # couldBeIncreased[s] is false if no unit could be added to the set s anymore
    couldBeIncreased = Vector{Bool}([])
    
    isInASet = Vector{Bool}(zeros(positionCount))
    
    # For each element in each initial set
    for setId in 1:length(sets)

        set = sets[setId]
        push!(hasCompatibilityWithSet, zeros(positionCount))
        push!(couldBeIncreased, true)

        for element in set
            # The element has compatibility with its set
            hasCompatibilityWithSet[setId][element] = true
            isInASet[element] = true
 
            for u in 1:positionCount
                if !areIncompatible[u, element]
                    hasCompatibilityWithSet[setId][u] = true
                end 
            end
        end
    end 
    
    # Id of the next set in which we will try to add a unit (or 0 if we try to create a new set)
    testedSetIdForAddition = 0

    isPossibleToCreateNewCluster = true

    # While all sets have not been tested to add units into them
    while testedSetIdForAddition <= length(sets)

        # If:
        # - we try to create a new cluster and we have not already tried to add all elements in a new cluster; or
        # - if we have not already tried to add all elements to the currently tested set and if their is at least two sets (if there is only one set, the inequality is satisfied by all fractional solutions)
        if testedSetIdForAddition == 0 && isPossibleToCreateNewCluster || couldBeIncreased[testedSetIdForAddition] && length(sets) > 1
            
            testedElementCount = 1
            elementsWeight = copy(initialElementsWeight.+1E-3)
            elementAdded = false

            # While we have not already tried to add all elements to the set and no unit has been added
            while !elementAdded && testedElementCount <= positionCount
                
                testedElement = sample(Weights(elementsWeight))

                # Set the weight to 0 to avoid picking the element again
                elementsWeight[testedElement] = 0

                if !isInASet[testedElement]
                    
                    # Test if the unit is incompatible with all the other sets
                    testedSetIdForIncompatibility = 1
                    isValidForAddition = true

                    while isValidForAddition && testedSetIdForIncompatibility <= length(sets)
                        if testedSetIdForIncompatibility != testedSetIdForAddition && hasCompatibilityWithSet[testedSetIdForIncompatibility][testedElement]
                            isValidForAddition = false
                        end

                        testedSetIdForIncompatibility += 1
                    end
                    # If the unit is not incompatible with all the other sets
                    if isValidForAddition
                        
                        # If the unit is added to a new cluster
                        if testedSetIdForAddition == 0
                            push!(sets, [testedElement])
                            push!(hasCompatibilityWithSet, zeros(positionCount))
                            push!(couldBeIncreased, 1, true)

                            # Temporarily set this to the new set id (it will be resetted to 0 before the next iteration)
                            testedSetIdForAddition = length(sets)
                        else
                            for element in sets[testedSetIdForAddition]
                              if areIncompatible[element, testedElement]
                                isValidForAddition = false
                                break
                              end
                            end
                            # Add the unit to the set
                            push!(sets[testedSetIdForAddition], testedElement)
                        end 

                        isInASet[testedElement] = true
                        
                        for u in 1:positionCount
                            if !areIncompatible[u, testedElement]
                                hasCompatibilityWithSet[testedSetIdForAddition][u] = true
                            end 
                        end

                        # We start again from the smallest set to add units
                        if isPossibleToCreateNewCluster
                            testedSetIdForAddition = 0
                        else
                            testedSetIdForAddition = 1
                        end 

                        elementAdded = true
                        
                    end # if isValidForAddition
                end # if !isInASet[testedElement]

                testedElementCount += 1
                
            end # while !elementAdded && testedElementCount <= positionCount

            if !elementAdded

                # If it was not possible to create a new cluster
                if testedSetIdForAddition == 0
                    isPossibleToCreateNewCluster = false                                
                else # If no unit can be added to an existing cluster
                    couldBeIncreased[testedSetIdForAddition] = false
                end 
                testedSetIdForAddition += 1
            end
        else
            testedSetIdForAddition += 1
        end # if testedSetIdForAddition == 0 && isPossibleToCreateNewCluster || couldBeIncreased[testedSetIdForAdd
    end # while testedSetIdForAddition < length(sets) || testedSetIdForAddition == length(sets) && length(sets) > 1 && length(sets[end]) == length(sets[end-1])

    if length(sets) > 1
        #println(length(sets), " set(s) in the incompatible sets inequalities: ", length.(sets))
    end 
    return sets
end 

"""
Build an IC inequality that corresponds to sets of elements
"""
function getICInequality(iFamily::InequalityFamily, model::Model, instance::Instance, sets::Vector{Vector{Int}}, dVariables::Dict{Symbol, Any})

    # If it is the first time 
    if !haskey(iFamily.additionalData, :previousClusters)
        previousClusters = Vector{Vector{Int}}()
        previousClustersVariables = Vector{VariableRef}()
        iFamily.additionalData[:previousClusters] = previousClusters
        iFamily.additionalData[:previousClustersVariables] = previousClustersVariables
    end

    # List of the clusters previously used in found constraints ordered lexicographically
    previousClusters = iFamily.additionalData[:previousClusters]

    # Variable which represent the covering of a cluster (ordered as variable previousClusters)
    previousClustersVariables = iFamily.additionalData[:previousClustersVariables]
    
    y, z, z2 = getVariables(dVariables)

    n = 0
    m = 0
    if z != nothing
        n = length(z)
    else
        m = size(z2, 1)
        n = size(z2, 2)
    end 
    
    if length(sets) > 1

        # Variable of each cluster
        constraintVariables = Vector{VariableRef}() 

        for setId in 1:length(sets)
            set = sets[setId]

            # If the set is of size 1, directly use variable z[u] or y[u]
            if length(set) == 1
                if set[1] <= n
                    if z != nothing 
                        push!(constraintVariables, z[set[1]])
                    else
                        for j in 1:m
                            push!(constraintVariables, z2[j, set[1]])
                        end 
                    end 
                else
                    push!(constraintVariables, y[set[1] - n])
                end 
            else # If the set has several elements
                sort!(set)
                insertIndex = searchsortedfirst(previousClusters, set)

                # If the set has already been used in a previous constraint, use the corresponding already created variable
                if insertIndex <= length(previousClusters) && previousClusters[insertIndex] == set
                    push!(constraintVariables, previousClustersVariables[insertIndex])
                    
                else # Otherwise, create a new variable for this set
                    newClusterVariable = @variable(model)
                    push!(constraintVariables, newClusterVariable)
                    
                    @constraint(model, newClusterVariable <= 1)

                    if z != nothing
                        @constraint(model, [u in set; u <= n], newClusterVariable >= z[u])
                    else
                        @constraint(model, [u in set; u <= n], newClusterVariable >= sum(z2[j, u] for j in 1:m))
                    end
                    
                    @constraint(model, [j in set; j >  n], newClusterVariable >= y[j - n]) 

                    # Update previousClusters and previousClustersVariables with this set and its variable
                    insert!(previousClusters, insertIndex, set)
                    insert!(previousClustersVariables, insertIndex, newClusterVariable)
                end 
            end 
        end
        return @build_constraint(sum(var for var in constraintVariables) <= 1)
    end
    println("Warning: no IC inequality added, only one set of elements provided")
    return nothing
end

"""
Add to a model an IC inequality that corresponds to sets of elements.
This is the weaker version of the inequality in which we do not add one variable to the model for each cluster.
The advantage is that it can then be used in a branch-and-cut (since it does not add variables to the model).
"""
function getWeakerICInequality(iFamily::InequalityFamily, model::Model, instance::Instance, sets::Vector{Vector{Int}}, dVariables::Dict{Symbol, Any})
    
    y, z, z2 = getVariables(dVariables)

    n = 0
    m = 0
    if z != nothing
        n = length(z)
    else
        m = size(z2, 1)
        n = size(z2, 2)
    end 
    
    if length(sets) > 1

        # Get the maximal length of a set...
        maxLength = 0

        ## ... and the variable associated to all elements in sets
        vars = Vector{VariableRef}()

        # For each set
        for set in sets

            # Update the maximal length of a set
            if length(set) > maxLength
                maxLength = length(set)
            end

            # For each element in the set add its corresponding variable(s) to array vars
            for element in set
                if element <= n
                    if z != nothing 
                        push!(vars, z[element])
                    else
                        for j in 1:m
                            push!(vars, z2[j, element])
                        end 
                    end 
                else
                    push!(vars, y[element - n])
                end
            end 
        end
        
        return @build_constraint(sum(var for var in vars) <= maxLength)
    end
    return nothing
end 

"""
Try to find an incompatible clusters inequality that separate a fractional solution
- dVariables: dictionary containing the variables used. It must contain one of the two next entries:
  - :z for single indexed z variables (z[u]=1 iff u is covered); or 
  - :z2 for doubled indexed z variables (z[j, u]=1 iff u is covered by site j); and it can additionally contain 
  - (optional) :y for the sites variables (y[j]=1 iff site j has an HAPS)

"""
function icSeparation(iFamily::InequalityFamily, instance::Instance, model::Model, dVariables::Dict{Symbol, Any}; cb_data=nothing)

    # Compute incompatibilities if it is the first time the separation is called
    if !haskey(iFamily.additionalData, :areIncompatible)

        positionCount = instance.n
        if haskey(dVariables, :y)
            positionCount = instance.n + instance.m
        end 

        iFamily.additionalData[:areIncompatible] = getIncompatibilities(instance, positionCount)
    end

    areIncompatible = iFamily.additionalData[:areIncompatible]

    # If there are incompatible elements
    if areIncompatible != nothing
        
        maxCutsAdded = min(iFamily.maxCutsAddedInOneIteration, iFamily.maxCutsAddedInNextIteration)
        maxCutsFound = 1.5 * maxCutsAdded
        
        y, z, z2 = getVariables(dVariables)
        positionCount = size(areIncompatible, 1)
        
        ## Set the weights to the value of the fractional solution
        elementsWeight = Vector{Float64}(zeros(positionCount))
        
        # If the formulation contains variables z_u
        if z != nothing
            elementsWeight[1:instance.n] = getValue.(z, cb_data) 
        else # If the formulation contains variables z_j,u
            for u in 1:instance.n
                elementsWeight[u] = maximum(getValue.(z2[:, u], cb_data))
            end
        end 
        
        if y != nothing
            elementsWeight[instance.n+1:end] = getValue.(y[:], cb_data)
        end 

        ## Get incompatible clusters
        cutsFound = Vector{Any}()

        # Value of the violation of each cut
        cutsViolation = Vector{Float64}()

        # The greater this weight for an element, the more likely it is to be picked in the next cut
        candidateWeight = copy(elementsWeight)
        remainingCandidates = count(x -> x>1E-4, elementsWeight)

        startingTime = time()
        elapsedTime = 0

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

        #@show round(iFamily.additionalData[:relaxationTime], digits = 4)
        #@show round(maxSeparationTimeWithoutCut, digits=2), round(maxSeparationTimeWithCut, digits=2)
        
        # While:
        # - the maximal number of cut found is not reached; and
        # - all the elements with non null variable value have not been used as seeds; and
        # - the separation time is lower than the relaxation time or no cut is found; and
        # - the separation time is lower than 5*(the relaxation time).
        while length(cutsFound) < maxCutsFound && remainingCandidates > 0 && (maxSeparationTimeWithCut == -1 || elapsedTime < maxSeparationTimeWithCut || length(cutsFound) == 0) && (maxSeparationTimeWithoutCut == -1 || elapsedTime < maxSeparationTimeWithoutCut) && (time_limit == -1 || elapsedTime <= time_limit)

            #@show elapsedTime
            seedElement = sample(Weights(candidateWeight))

            sets = generateIC(instance, areIncompatible, initialElementsWeight=elementsWeight, initialSets = Vector{Vector{Int}}([[seedElement]]))

            cutFound, violation = isViolatedByIC(model, dVariables, sets, cb_data = cb_data)

            if cutFound                
                push!(cutsFound, sets)
                push!(cutsViolation, violation)
                #@show seedElement, remainingCandidates, sets#, elementsWeight 

                for set in sets
                    for element in set
                        if candidateWeight[element] > 0
                            remainingCandidates -= 1
                        end 
                        candidateWeight[element] = 0
                    end
                end
            else
                candidateWeight[seedElement] = 0
                remainingCandidates -= 1
            end
            elapsedTime = time() - startingTime
        end

        #@show length(cutsFound), round(elapsedTime, digits=2)
        #readline()

        
        cutsFound = cutsFound[sortperm(cutsViolation, rev = true)]
        return cutsFound[1:min(end, maxCutsAdded)]
    else
        return Vector{Any}([])
    end
end 

"""
Test if a fractional solution is violated by and incompatible clusters inequalities
"""
function isViolatedByIC(model::Model, dVariables::Dict{Symbol, Any}, sets::Vector{Vector{Int}}; cb_data=nothing)

    y, z, z2 = getVariables(dVariables)
    singleIndexedZVariables = z != nothing
    lhs = 0.0
    n = 0

    if singleIndexedZVariables
        n = length(z)
    else
        n = size(z2, 2)
    end 
    
    for set in sets

        maxSetValue = 0.0 
        for element in set

            if element <= n
                maxZValue = nothing
                if singleIndexedZVariables
                    maxZValue = getValue(z[element], cb_data)
                else
                    maxZValue = maximum(getValue.(z2[:, element], cb_data))
                end
                
                if maxZValue > maxSetValue
                    maxSetValue = maxZValue
                end
            else
                if getValue(y[element - n], cb_data) > maxSetValue
                    maxSetValue = getValue(y[element - n], cb_data)
                end
            end
        end

        lhs += maxSetValue
    end

    violation = lhs - 1

    return violation > 1E-4, violation
end

function getIncompatibilities(instance, positionCount)

    maximalClientsDistance      = (instance.p - 1) * instance.rCom + 2 * sqrt(instance.rCouv^2 - instance.L^2) 
    maximalSitesDistance        = (instance.p - 1) * instance.rCom
    # Wrong bound (assumes that the distance between a site and a client does not take into account their differeng heights)
    # maximalClientToSiteDistance = (instance.p - 1) * instance.rCom + sqrt(instance.rCouv^2 - instance.L^2)

    maximalClientToSiteDistance = sqrt(((instance.p - 1) * instance.rCom + sqrt(instance.rCouv^2 - instance.L^2))^2+L^2)
    
    xMin = minimum(instance.clientPositions[:, 1]) 
    xMax = maximum(instance.clientPositions[:, 1])
    yMin = minimum(instance.clientPositions[:, 2])
    yMax = maximum(instance.clientPositions[:, 2])

    compatibilityFound = false

    if maximalClientsDistance < sqrt((xMax-xMin)^2+(yMax-yMin)^2)

        areIncompatible = Matrix{Bool}(zeros(positionCount, positionCount))
        
        for u1 in 1:instance.n
            for u2 in u1+1:instance.n
                dist = sqrt((instance.clientPositions[u1, 1]-instance.clientPositions[u2, 1])^2 + (instance.clientPositions[u1, 2]-instance.clientPositions[u2, 2])^2)
                if dist > maximalClientsDistance
                    compatibilityFound = true
                    areIncompatible[u1, u2] = true
                    areIncompatible[u2, u1] = true
                end 
            end

            if positionCount > instance.n
                for j in 1:instance.m
                    if instance.d[u1, j] > maximalClientToSiteDistance
                        compatibilityFound = true 
                        areIncompatible[u1, instance.n + j] = true
                        areIncompatible[instance.n + j, u1] = true
                    end 
                end
            end 
        end

        if positionCount > instance.n
            for j1 in 1:instance.m
                for j2 in j1+1:instance.m
                    if instance.d_com[j1, j2] > maximalSitesDistance
                        compatibilityFound = true 
                        areIncompatible[instance.n + j1, instance.n + j2] = true
                        areIncompatible[instance.n + j2, instance.n + j1] = true
                    end 
                end
            end
        end 
    else
        println("No incompatibilities: the maximal distance between two covered clients is larger than the distance between any pair of clients. Non incompatible set inequalities will be generated.")
        areIncompatible = nothing
    end

    if !compatibilityFound
        println("No incompatibilities found. Non incompatible set inequalities will be generated.")
        areIncompatible = nothing
    end

    return areIncompatible
end 

function cleanPreviousClusters(iFamily::InequalityFamily)
    iFamily.additionalData[:previousClusters] = Vector{Vector{Int}}()
    iFamily.additionalData[:previousClustersVariables] = Vector{VariableRef}()
end 
