using Graphs
#include("../distance_computation_3D.jl")

"""
An instance of the localisation problem
"""
mutable struct Instance

    # Number of clients
    n::Int

    # Number of candidate sites
    m::Int

    # Maximal radius at which a drone can reach a client
    rCouv::Int

    # Maximal radius at which two drones can communicate
    rCom::Int

    # Height of the drones
    L::Int

    # Number of drones used
    p::Int

    # Distance between the clients and the candidate sites
    # d[i, j] = distance between client i and site j
    d::Array{Float64, 2}

    # Distances between the candidate sites
    d_com::Array{Float64, 2}

    # Position of the clients
    clientPositions::Array{Float64, 2}

    # Positions of the sites
    # sitesPositions[j, 1]: x-axis coordinate of position j
    # sitesPositions[j, 2]: y-axis coordinate of position j
    sitesPositions::Array{Float64, 2}

    # areClientsClose[i, j] is true if client i and j can be covered by the same HAPS
    areClientsClose::Matrix{Bool}
    
    function Instance()
        return new()
    end   
end 
  
"""
Instance constructor which compute the distances
"""
function Instance(n::Int, m::Int, rCouv::Int, rCom::Int, L::Int, p::Int, clientPositions::Array{Float64, 2}, sitesPositions::Array{Float64, 2}; clientHeight::Vector{Float64}=Vector{Float64}([]))  
    clientPositions = round.(clientPositions, digits = 2)
    sitesPositions = round.(sitesPositions, digits = 2)
    # Compute the distance between clients and sites
    d = Matrix{Float64}([0 for i in 1:n, j in 1:m])
    #println("Warning: client height ignored (not yet implemented in quadratic formulations)") 
    for i in 1:n, j in 1:m
#        if clientHeight == []
            d[i,j]=sqrt( (clientPositions[i, 1]-sitesPositions[j, 1])^2 + (clientPositions[i, 2]-sitesPositions[j, 2])^2 + L^2 )
#        else
#            d[i,j]=round(sqrt( (clientPositions[i, 1]-sitesPositions[j, 1])^2 + (clientPositions[i, 2]-sitesPositions[j, 2])^2 + (L-clientHeight[i])^2))
#        end 
    end

    # Compute the distance between sites
    d_com = Matrix{Float64}([0 for j in 1:m, j1 in 1:m])
    for j in 1:m, j1 in 1:m
        d_com[j,j1]=sqrt( (sitesPositions[j, 1]-sitesPositions[j1, 1])^2 + (sitesPositions[j, 2]-sitesPositions[j1, 2])^2 )
    end

    return Instance(n, m, rCouv, rCom, L, p, d, d_com, clientPositions, sitesPositions)
end

"""
Instance constructor which does not compute the distances
"""
function Instance(n::Int, m::Int, rCouv::Int, rCom::Int, L::Int, p::Int, d::Array{Float64, 2}, d_com::Array{Float64, 2}, clientPositions::Array{Float64, 2}, sitesPositions::Array{Float64, 2})   
  
    this = Instance()  
    this.n = n
    this.m = m
    this.rCouv = rCouv
    this.rCom = rCom
    this.L = L
    this.p = p
    this.clientPositions = clientPositions
    this.sitesPositions = sitesPositions
    this.d = d
    this.d_com = d_com
    this.areClientsClose = Matrix{Bool}(zeros(0, 0))
  
    return this   
end

"""
Instance constructor from the path of a file
"""
function Instance(path::String; recomputeDistances::Bool=false)
    global clientHeight = []

    d = nothing
    d_com = nothing

    include(path)
    
    if d == nothing || d_com == nothing
        recomputeDistances = true
    end 
    
    if clientHeight == [] && !recomputeDistances

        return Instance(n, m, rCouv, rCom, L, p, d, d_com, clientPositions, sitesPositions)
    else
        if clientHeight != []
            instance =  Instance(n, m, rCouv, rCom, L, p, clientPositions, sitesPositions, clientHeight = clientHeight)
        else
            instance =  Instance(n, m, rCouv, rCom, L, p, clientPositions, sitesPositions)
        end 

        if recomputeDistances
            open(path, "w") do fout
                println(fout, instance)
            end
        end
        return instance
    end 
end 

"""
How to print the structure
(used to write instance files)
"""
function Base.show(io::IO, instance::Instance)    
    println(io, "n = ", instance.n)
    println(io, "m = ", instance.m)
    println(io, "rCouv = ", instance.rCouv)
    println(io, "rCom = ", instance.rCom)
    println(io, "L = ", instance.L)
    println(io, "p = ", instance.p)
    println(io, "d = ", instance.d)
    println(io, "d_com = ", instance.d_com)
    println(io, "clientPositions = ", instance.clientPositions)
    println(io, "sitesPositions = ", instance.sitesPositions)
end   

"""
Generate an instance and save it in a text file

The position of the clients and the sites are randomly generated in [0, 100]^2.

Input:
- n, m, rCOuv, rCom, L, p: parameters of the instance
- (optional) outputFile: path in which the instance is saved ("" if no file is saved)
"""
function generateInstance(n::Int, m::Int, rCouv::Int, rCom::Int, L::Int, p::Int; outputFile::String="")

    # Clients coordinates
    clientPositions=[100*rand() for i in 1:n, k in 1:2] 
    
    # Sites coordinates
    sitePositions=[100*rand() for j in 1:m, k in 1:2]
    instance = Instance(n, m, rCouv, rCom, L, p, clientPositions, sitePositions) 

    if outputFile != ""
        open(outputFile, "w") do fout
            println(fout, instance)
        end
    end

    return instance
end

"""
Generate a set of instance with the same unit position but in which the position of the sites are grid of different length
"""
function generateGridInstances(n::Int, fieldLength::Int, gridLengths::Vector{Int}, rCouv::Int, rCom::Int, L::Int, p::Int, outputPath::String=""; altitudeMap=nothing, saveClientSiteDistances::Bool=true)

    if !isdir(outputPath)
        mkpath(outputPath)
    end
    
    # Clients coordinates
    clientPositions = [fieldLength*rand() for i in 1:n, k in 1:2]

    clientHeight = Vector{Float64}(zeros(n))
    
    # For each grid length considered
    for gridLength in gridLengths
        
        outputFile = outputPath * "n" * string(n) * "_gridLength" * string(gridLength) * "_rCom" * string(rCom) * ".txt"
        
        # Sites coordinates        
        sitePositions = generateGrid(fieldLength, fieldLength, gridLength)

        instance = Instance(n, size(sitePositions, 1), rCouv, rCom, L, p, clientPositions, sitePositions)

        if altitudeMap != nothing
            instance.d = get_clients_sites_distances(instance.clientPositions, instance.sitesPositions, fieldLength, altitudeMap, L)

            if !saveClientSiteDistances
                for i in 1:n
                    clientHeight[i] = get_point_altitude(instance.clientPositions[i, 1], instance.clientPositions[i, 2], fieldLength, altitudeMap)
                end 
            end 
        end 

        if saveClientSiteDistances
            open(outputFile, "w") do fout
                println(fout, instance)
            end
        else
            open(outputFile, "w") do fout
                println(fout, "n = ", instance.n)
                println(fout, "m = ", instance.m)
                println(fout, "rCouv = ", instance.rCouv)
                println(fout, "rCom = ", instance.rCom)
                println(fout, "L = ", instance.L)
                println(fout, "p = ", instance.p)
                println(fout, "d = ", [])
                println(fout, "d_com = ", [])
                println(fout, "clientPositions = ", instance.clientPositions)
                println(fout, "sitesPositions = ", instance.sitesPositions)
                println(fout, "clientHeight = ", clientHeight)
            end
        end 
    end 
end 

"""
Generate a grid of positions.

Input:
- fieldXSize: horizontal length of the field
- fieldYSize: vertical length of the field
- gridLength: distance between two adjacent positions in the grid

Output:
- gridPositions::Matrix{Float64} such that gridPositions[i, j] is the jth coordinate of the ith position
"""
function generateGrid(fieldXSize::Int, fieldYSize::Int, gridLength::Int)
    
    nextX = 0
    nextY = 0

    gridPositions = Matrix{Float64}(zeros(0, 2))

    while nextY <= fieldYSize
        gridPositions = vcat(gridPositions, [nextX nextY])

        nextX += gridLength
        
        if nextX > fieldXSize
            nextX = 0
            nextY += gridLength
        end 
    end

    return gridPositions
end 


function generateDataset()

    for n in 10:10:100
        for p in [3, 6, 9]
            rCom = 50

            if p == 6
                rCom /= 1.5
            elseif p == 9
                rCom /= 2
            end
            rCom = round(Int, rCom)

            outputFile = "../../data/expe1/n" * string(n) * "_m" * string(n) * "_p" * string(p) * "_rCom" * string(rCom)  * ".txt"

            if !isfile(outputFile)
                generateInstance(n, n, rCom, rCom, 20, p, outputFile = outputFile)
            end 
        end
    end 
end 



function generateCAIDDataset()

    gridLengths = [125, 250, 500, 1000]
    rCom = 1000
    rCouv = 1000
    L = 500
    fieldLength = 5000
    
    outputPath = "../data/expeCAID1/set2/" 
    altitudeMap = getBrunoAltitudeMap()
    for n in 25:25:100
        generateGridInstances(n, fieldLength, gridLengths, rCom, rCouv, L, 1, outputPath, altitudeMap = altitudeMap) 
    end 
end

function generateCAIDDataset4()

    gridLengths = [100, 150, 200, 300, 375, 400, 500, 600, 750, 1000, 1200, 1500]
    rCom = 1500
    rCouv = 1500
    L = 1200
    fieldLength = 6000
    
    outputPath = "../data/expeCAID1/set4/" 
    altitudeMap = getBrunoAltitudeMap()
    for n in [100, 1000]
        generateGridInstances(n, fieldLength, gridLengths, rCom, rCouv, L, 1, outputPath, altitudeMap = altitudeMap) 
    end 
end 


function generateCAIDDataset5()

    gridLengths = [100, 150, 200, 300, 375, 400, 500, 600, 750, 1000, 1200, 1500]
    rCom = 1500
    rCouv = 1500
    L = 1200
    fieldLength = 6000
    saveClientSiteDistances = false
    
    outputPath = "../data/expeCAID1/set5/" 
    altitudeMap = getBrunoAltitudeMap()
    for n in [100, 1000]
        generateGridInstances(n, fieldLength, gridLengths, rCom, rCouv, L, 1, outputPath, altitudeMap = altitudeMap, saveClientSiteDistances = saveClientSiteDistances) 
    end 
end 


function generateCAIDDataset7()

    gridLengths = [200, 300, 375, 400, 500, 600, 750, 1000, 1200, 1500]
    rCom = 1500
    rCouv = 1500
    L = 1200
    fieldLength = 6000
    saveClientSiteDistances = false
    
    outputPath = "../data/expeCAID1/set7/" 
    altitudeMap = getBrunoAltitudeMapMax120()
    for n in [100, 1000]
        generateGridInstances(n, fieldLength, gridLengths, rCom, rCouv, L, 1, outputPath, altitudeMap = altitudeMap, saveClientSiteDistances = saveClientSiteDistances) 
    end 
end

function generateCAIDDataset8()

    gridLengths = [200, 300, 375, 400, 500, 600, 750, 1000, 1200, 1500]
    rCom = 1500
    rCouv = 1500
    L = 1200
    fieldLength = 6000
    saveClientSiteDistances = false
    
    outputPath = "../data/expeCAID1/set7/" 
    altitudeMap = getBrunoAltitudeMapMax120()
    for n in [100, 1000]
        generateGridInstances(n, fieldLength, gridLengths, rCom, rCouv, L, 1, outputPath, altitudeMap = altitudeMap, saveClientSiteDistances = saveClientSiteDistances) 
    end 
end 

function generateLargeDataset()

    gridLengths = [50, 100, 200]
    rCom = 1500
    rCouv = 1500
    L = 1200
    fieldLength = 6000
    saveClientSiteDistances = false
    
    outputPath = "../data/2022-12-large_instances/" 
    for n in [2000, 5000, 10000]
        generateGridInstances(n, fieldLength, gridLengths, rCom, rCouv, L, 1, outputPath, saveClientSiteDistances = saveClientSiteDistances) 
    end 
end 

function randomizeSitePositions(instancePath::String)
    
    fieldLength = 6000
    instance = Instance(instancePath)

    include(instancePath)

    sitesPositions = [fieldLength*rand() for i in 1:instance.m, k in 1:2]

    outputPath = replace(instancePath, ".txt" => "_randomized.txt")
    
    open(outputPath, "w") do fout
        println(fout, "n = ", instance.n)
        println(fout, "m = ", instance.m)
        println(fout, "rCouv = ", instance.rCouv)
        println(fout, "rCom = ", instance.rCom)
        println(fout, "L = ", instance.L)
        println(fout, "p = ", instance.p)
        println(fout, "d = ", [])
        println(fout, "d_com = ", [])
        println(fout, "clientPositions = ", instance.clientPositions)
        println(fout, "sitesPositions = ", sitesPositions)
        println(fout, "clientHeight = ", clientHeight)
    end
    
end 

function addAltitudeToInstance(path::String)
    instance = Instance(path)

    altitudeMap = getBrunoAltitudeMap()
    instance.d = get_clients_sites_distances(instance.clientPositions, instance.sitesPositions, 6000, altitudeMap, 1200) 
    outputPath = "../data/expeCAID3/" 

    mkpath(outputPath)
    open(outputPath * basename(path), "w") do fout
        println(fout, instance)
    end
end 

function computeCloseClients!(instance)

    println("Computing close clients...")
    startingTime = time()
    instance.areClientsClose = Matrix{Bool}(zeros(instance.n, instance.n))

    closeThreshold = 2 * sqrt(instance.rCouv^2 - instance.L^2)###consider alttitude
    for u1 in 1:instance.n
        for u2 in u1+1:instance.n###no repeat
            dist = sqrt((instance.clientPositions[u1, 1] - instance.clientPositions[u2, 1])^2+(instance.clientPositions[u1, 2] - instance.clientPositions[u2, 2])^2)
            if dist <= closeThreshold
                instance.areClientsClose[u1, u2] = true
                instance.areClientsClose[u2, u1] = true
            end 
        end
    end

    println("Computing close clients done in ", round(Int, time() - startingTime),  "s")
end

function getMaximalCliques!(instance)
    
    if size(instance.areClientsClose, 1) == 0
        computeCloseClients!(instance)
    end

    g = SimpleGraph()
    add_vertices!(g, instance.n)

    for u1 in 1:instance.n
        for u2 in u1+1:instance.n
            if !instance.areClientsClose[u1, u2]
                add_edge!(g, u1, u2)
            end
        end
    end

    println("Computing maximal cliques...")
    cliques = maximal_cliques(g)
    println("Computing maximal cliques done in ", round(Int, time() - startingTime),  "s")
    
    return cliques
    
end

function getCliquesHeuristically!(instance::Instance)

    if size(instance.areClientsClose, 1) == 0
        computeCloseClients!(instance)
    end

    degree = Vector{Int}(zeros(instance.n))
    
    for u in 1:instance.n
        degree[u] = sum(instance.areClientsClose[u, :] .== false)
    end

    # isInAClique[u] is true if client u is in one of the found cliques
    isInAClique = Vector{Bool}(zeros(instance.n))

    sortedDegrees = sortperm(degree, rev=true)

    clique = Vector{Int}([])

    isValidCandidate = Vector{Bool}(ones(instance.n))

    ## 1 - For each client in decreasing order of their degrees, try to add it to a first clique
    for clientDegreeId in 1:instance.n

        clientId = sortedDegrees[clientDegreeId]

        if isValidCandidate[clientId]
            push!(clique, clientId)

            for clientDegreeId2 in clientDegreeId+1:instance.n
                clientId2 = sortedDegrees[clientDegreeId2]
                if instance.areClientsClose[clientId, clientId2]
                    isValidCandidate[clientId2] = false
                end 
            end 
        end 
    end

    cliques = Vector{Vector{Int}}([])

    if length(clique) > instance.p
        push!(cliques, clique)
        for client in clique
            isInAClique[client] = true
        end 
    end 

    ## 2 - For each client which is not in a clique, try to create a clique which includes it by randomly choosing clients
    clientDegreeId = 1

    # For each client...
    while clientDegreeId <= instance.n

        clientId = sortedDegrees[clientDegreeId]

        # If the client has enough neighbor to lead to a clique of size p+1
        if degree[clientId] > instance.p

            # If the client is not already in a clique
            if !isInAClique[clientId] && degree[clientId] > instance.p

                newClique = Vector{Int}([clientId])

                # Randomly pick clients which are connected to all the clients in the new clique
                # (clients which have a high degree will have a greater chance of being picked)
                clientWeight = copy(degree)
                
                candidateClientsCount = instance.n
                
                # All clients wich are close to clientId or do not have enough neighbors have a weight of 0
                for clientId2 in 1:instance.n
                    if instance.areClientsClose[clientId, clientId2] || clientId == clientId2 || degree[clientId2] < instance.p + 1
                        clientWeight[clientId2] = 0
                        candidateClientsCount -= 1
                    end
                end
                
                while candidateClientsCount > 0

                    # Randomly pick a client with a non-null weight and add it to the clique
                    nextClient = sample(Weights(clientWeight))
                    push!(newClique, nextClient)
                    
                    # Set the weight of each client which is close to the new client to 0
                    # idem if the client does not have enough neighbors to reach a clique of size >= p+1
                    for clientId2 in 1:instance.n
                        if instance.areClientsClose[nextClient, clientId2] || clientId2 == nextClient || degree[clientId2] < instance.p + 1 - length(newClique)
                            clientWeight[clientId2] = 0
                            candidateClientsCount -= 1
                        end
                    end 
                end

                if length(newClique) > instance.p
                    push!(cliques, newClique)

                    for clientId2 in newClique
                        isInAClique[clientId2] = true
                    end 
                end 
            end
        else # If there is no more clients with enough neighbors
            clientDegreeId = instance.n
        end
        
        clientDegreeId += 1 
    end 

    return cliques
end


function generateInstanceWithGap()

    n = 10
    m = 8
    rCouv = 30
    rCom = 15
    L = 20
    p = 3

    res = nothing

    param = Dict{String, Any}()
    param["redundancyRatio"]    = 0.5
    param["isRelaxation"]       = false
    param["threads"]            = -1
    param["addMaxCliques"]      = true
    param["tightenRelaxation"]  = false
    param["time_limit"]         = -1
    param["warmStart"]          = true
    param["maxCutsCount"]       = 10
    param["connexity"]          = "MTZ"
    param["instancePath"]       = "instanceWithGap.txt"
    param["preprocessing"]      = false
    param["p"]                  = 3
    instance = nothing

    while res == nothing || res["rootRelaxation"] - res["objective"] < 0.4
        instance = generateInstance(n, m, rCouv, rCom, L, p; outputFile="instanceWithGap.txt")
        res = aeSolve(param)
    end
end 




