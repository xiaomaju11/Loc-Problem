
mutable struct InequalityFamily

    # Each element of this vector represents a cut added to the model
    # (in such a way that the cut associated to an element can be reconstructed from the element)
    addedCuts::Vector{Any}
    
    # Function that cuts a fractional point using an inequality of the current family
    # The violated inequalities are added to the model in this function.
    #
    # Signature: fSeparate(iFamily::InequalityFamily, instance::Instance, model::Model, dVariables::Dict{Symbol, Any}; cb_data=nothing)
    #
    # Inputs
    # 1. the inequality family
    # 2. instance: the instance
    # 3. model: the model which current solution is the fractional point to cut
    # 4. dVariables: the variables of the model identified by a symbol
    # 5. cb_data: data of the callback or nothing if we are in a cutting plane
    #
    # Output
    # results::Vector{Any}: vector of violated cuts (as represented in this.addedCuts) added during the separation
    fSeparation::Function

    # Function that generate a cut of this family
    #
    # Signature: fGetCut(inequalityFamily::InequalityFamily, model::Model, instance::Instance, separationResult::Any, dVariables::Dict{Symbol, Any})
    #
    # Inputs
    # 1. inequalityFamily: the current inequality family
    # 2. model: the model in which the cut must be added
    # 3. instance: the instance currently solved
    # 4. separationResult: one cut (as represented in this.addedCuts)
    # 5. dVariables: the variables of the model
    fGetCut::Function

    # Maximal number of cuts added at each call of the separation procedure
    maxCutsAddedInOneIteration::Int

    # Maximal number of cuts added in the next separation
    maxCutsAddedInNextIteration::Int

    # Any data required to separate the inequality of this family
    # (e.g., incompatible clusters requires to know which units are not compatible which should not be recomputed at each execution of the separation procedure)
    # It also contains an entry :time_limit equal to -1 if the time is not limited
    additionalData::Dict{Symbol, Any}

    # Function 
    
    # Function called at the beginning and the end of the cutting plane (do nothing by default)
    # It can for example remove additional data that would not be relevant from one execution of the cutting plane to another (e.g., if the instance change)
    #
    # Signature: fClean(iFamily::InequalityFamily)
    fClean::Function

    function InequalityFamily()
        return new()
    end
end

function InequalityFamily(fSeparation, fGetCut; maxCutsAddedInOneIteration::Int64=typemax(Int64), fClean::Function=cleanNothing)

    this = InequalityFamily()

    this.fSeparation = fSeparation
    this.fGetCut = fGetCut
    this.maxCutsAddedInOneIteration = maxCutsAddedInOneIteration
    this.maxCutsAddedInNextIteration = typemax(Int64)
    this.addedCuts = Vector{Any}()
    this.additionalData = Dict{Symbol, Any}()
    this.fClean = fClean

    return this
end 

include("../incompatibleClustersInequalities.jl")
include("../cliqueInequalities.jl")

"""
Apply a cutting plane algorithm.

Input
- model: the continuous model to solve
- inequalityFamilies: the family of inequalities that will be separated in the cutting plane algorithm. Each family is a Tuple(Function, Function) in which:
  - the first function 
- (optional) maxCutsAdded: maximal number of cuts added in the cutting plane (-1 if it is not limited)
"""
function cuttingPlane(instance::Instance, model::Model, inequalityFamilies::Vector{InequalityFamily}, dVariables::Dict{Symbol, Any}; maxCutsAdded::Int64=-1, time_limit::Int=-1)

    cpStartingTime = time()
    totalSeparationTime = 0
    
    inequalityFound = true

    for family in inequalityFamilies
        family.fClean(family)
        family.addedCuts = Vector{Any}()
    end

    cutsAdded = 0
    firstRootRelaxationValue = nothing
    lastRootRelaxationValue = nothing

    firstRelaxationSolved = false
    relaxationSolved = true

    elapsedTime = time() - cpStartingTime
    
    while inequalityFound && cutsAdded < maxCutsAdded && relaxationSolved && (time_limit == -1 || elapsedTime <= time_limit)

        # Solve the relaxation
        relaxationStartingTime = time()

        if time_limit != -1
            set_time_limit_sec(model, max(5, time_limit - elapsedTime))
        end 
        optimize!(model)
        relaxationTime = time() - relaxationStartingTime
        solveRelaxation = termination_status(model) == MOI.NUMERICAL_ERROR

        # Give its time to all family of inequalities if it is the first relaxation
        if !firstRelaxationSolved
            for iFamily in inequalityFamilies
                iFamily.additionalData[:relaxationTime] = relaxationTime
            end 
            firstRelaxationSolved = true
        end 

        # Test if the relaxation is solved optimally (which may not be the case due to the time limit)
        relaxationSolved = termination_status(model) == MOI.OPTIMAL

        if relaxationSolved
            #displaySolution(dVariables, instance, isRelaxation=true)
            #readline()

            lastRootRelaxationValue = JuMP.objective_value(model)
            if firstRootRelaxationValue == nothing
                firstRootRelaxationValue = lastRootRelaxationValue
            end 
            
            println("--Relaxation value: ", round(lastRootRelaxationValue, digits = 2))
            #displaySolution(dVariables, instance, isRelaxation=true)
            
            inequalityFound = false
            familyKeyId = 1

            startingSeparationIterationTime = time()
            while familyKeyId <= length(inequalityFamilies) && !inequalityFound && (time_limit == -1 || elapsedTime <= time_limit)

                iFamily = inequalityFamilies[familyKeyId]
                if time_limit == -1
                    iFamily.additionalData[:time_limit] = -1
                else 
                    iFamily.additionalData[:time_limit] = time_limit - (time() - startingSeparationIterationTime)
                end 

                # Try to find violated inequalities and add them
                # (fSeparate must return an array in which each element corresponds to an added violated cut and an empty array if there are none)
                cutsDescription = iFamily.fSeparation(iFamily, instance, model, dVariables)

                for cut in cutsDescription
                    cstr = iFamily.fGetCut(iFamily, model, instance, cut, dVariables)
                    
                    if cstr != nothing
                        add_constraint(model, cstr)
                    end 
                end 
                append!(iFamily.addedCuts, cutsDescription)

                cutsAdded += length(cutsDescription)
                inequalityFound = length(cutsDescription) > 0

                #@show familyKeyId, cutsDescription
                familyKeyId += 1
                elapsedTime = time() - cpStartingTime
            end
            totalSeparationTime += time() - startingSeparationIterationTime

            #println("Elapsed cp time: ", round(Int, time() - cpStartingTime), "s")
            #println("Elapsed separation time: ", round(Int, totalSeparationTime), "s")
        end 
        elapsedTime = time() - cpStartingTime 
    end

    if inequalityFound
        optimize!(model)

        if termination_status(model) == MOI.OPTIMAL
            lastRootRelaxationValue = JuMP.objective_value(model)
        end 
    end 

    for family in inequalityFamilies
        family.fClean(family)
    end

    cpTime = time() - cpStartingTime

    #readline()
    return cpTime, totalSeparationTime, firstRootRelaxationValue, lastRootRelaxationValue, cutsAdded
    
end
function addCallback(instance::Instance, model::Model, inequalityFamilies::Vector{InequalityFamily}, dVariables::Dict{Symbol, Any})

    savedNodeId = -1
    savedDualBound = -1
    iterationsSinceLastSave = 0

    function cbFunction(cb_data::CPLEX.CallbackContext, context_id::Clong)

        if context_id != CPX_CALLBACKCONTEXT_RELAXATION
            return
        end
        println("Callback called:", length(inequalityFamilies[1].addedCuts), " cuts added")
        CPLEX.load_callback_variable_primal(cb_data, context_id)

        # Get the objective value of the fractional solution
        currentObjectiveValue = -1

        if haskey(dVariables, :z)
            currentObjectiveValue = sum(callback_value.(cb_data, dVariables[:z]))
        else
            currentObjectiveValue = sum(callback_value.(cb_data, dVariables[:z2]))
        end 


        # Get the id of the current B&B node
        valueP = Ref{Int64}()
        temp = CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODEUID, valueP)
        currentNodeId = valueP[]

        # If this is the first time the callback is called on this node
        if currentNodeId != savedNodeId
            savedNodeId = currentNodeId
            savedDualBound = currentObjectiveValue
            iterationsSinceLastSave = 1
        else
            iterationsSinceLastSave += 1

            # If the dual bound has not been sufficiently improved in the last 5 calls of the callback, stop
            if iterationsSinceLastSave > 5
                relaxationImprovement = abs(savedDualBound - currentObjectiveValue) / (savedDualBound+1E-4)

                if relaxationImprovement < 0.01
                    return 
                else
                    savedDualBound = currentObjectiveValue
                    iterationsSinceLastSave = 1
                end 
            end 
        end 
        # For each family of inequalities
        for iFamily in inequalityFamilies
            
            # Find cut(s)
            cutsDescription=iFamily.fSeparation(iFamily, instance, model, dVariables; cb_data=cb_data)
            #println("cutsDescription: ", cutsDescription)
            # And add them
            for cut in cutsDescription
                #MOI.set(model, MOI.UserCutCallback(), cb_data -> begin
                 cstr = iFamily.fGetCut(iFamily, model, instance, cut, dVariables)
                 #println("Generated constraint: ", cstr)
                 MOI.submit(model, MOI.UserCut(cb_data), cstr)
            end  
            append!(iFamily.addedCuts, cutsDescription)         
        end 
    end

    
    MOI.set(model, MOI.NumberOfThreads(), 1)
    # Add the callback to the model
    MOI.set(model, CPLEX.CallbackFunction(), cbFunction)
    
end


# Used in callbacks
function getWeakerICIFamily(maxCutsAddedInOneIteration::Int)
    return InequalityFamily(icSeparation, getWeakerICInequality; maxCutsAddedInOneIteration=maxCutsAddedInOneIteration, fClean=cleanPreviousClusters) 
end 

# Used in cutting plane
function getICIFamily(maxCutsAddedInOneIteration::Int)
    return InequalityFamily(icSeparation, getICInequality;maxCutsAddedInOneIteration=maxCutsAddedInOneIteration, fClean=cleanPreviousClusters) 
end 

function getCliqueIFamily(maxCutsAddedInOneIteration::Int)
    return InequalityFamily(cliqueSeparation, getCliqueInequality;maxCutsAddedInOneIteration=maxCutsAddedInOneIteration) 
end 

function cleanNothing(iFamily::InequalityFamily)
end 

function getValue(var::VariableRef, cb_data)
    if cb_data == nothing
        return JuMP.value(var)
    else
        callback_value(cb_data, var)
    end 
end 
