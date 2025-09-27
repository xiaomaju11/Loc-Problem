include("inequalityFamily.jl")

mutable struct ResolutionParam

    connexity::Symbol
    
    isRelaxation::Bool
    solveRelaxation::Bool
    displaySolution::Bool
    addMaxCliques::Bool
    tightenRelaxation::Bool
    preprocessing::Bool
    warmStart::Bool
    ciIneq::Bool
    addMaximalCliques::Bool
    useGridInitialization::Bool
    addEigenvalue::Bool

    verbose::Int64
    threads::Int64
    maxCuts::Int64
    maxCutsCount::Int64
    cpMaxCutsCount::Int64
    time_limit::Int64
    p::Int64
    
    redundancyRatio::Float64

    instancePath::String

    # Family of inequalities that are separated at the root node
    cpIneqFamilies::Vector{InequalityFamily}

    # Family of inequalities that can be separated at each node
    bcIneqFamilies::Vector{InequalityFamily}
    
    function ResolutionParam()
        return new()
    end
end

function ResolutionParam(connexity::String, instancePath::String; p::Int=2, time_limit::Int=-1)

    this = ResolutionParam()

    this.instancePath = instancePath
    this.connexity = Symbol(connexity)
    this.p = p
    this.time_limit = time_limit

    # Default values
    this.cpIneqFamilies = Vector{InequalityFamily}()
    this.bcIneqFamilies = Vector{InequalityFamily}()
    this.solveRelaxation = true
    this.isRelaxation = false
    this.threads = -1
    this.maxCuts = 50
    this.maxCutsCount = 0
    this.ciIneq = false
    this.maxCutsCount = 0 
    this.cpMaxCutsCount = 0 
    this.maxCuts = 50 
    this.preprocessing = false
    this.tightenRelaxation = false
    this.redundancyRatio = Float64(0.0)
    this.addMaxCliques = false
    this.displaySolution = false
    this.warmStart = true
    this.verbose = 1
    this.addMaximalCliques = false
    this.useGridInitialization = false
    this.addEigenvalue = false
    return this
end 
