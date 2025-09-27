using Dates

"""
Create a tex file which represents the solution included in a dictionary returned by a resolution method
"""
function texDocument(results::Dict{String, Any}; path::String="../results/output.tex")

    if haskey(results, "instancePath") && haskey(results, "isClientCovered") && haskey(results, "isSiteOpened")
        texDocument(Instance(results["instancePath"]), clientCoverage = results["isClientCovered"], siteOpening = results["isSiteOpened"], path = path)
        println("Tex file create in file ", path)
    elseif haskey(results, "instancePath") && haskey(results, "isClientCovered") && haskey(results, "xHAPS") && haskey(results, "yHAPS")
        texDocument(Instance(results["instancePath"]), clientCoverage = results["isClientCovered"], xHAPS = results["xHAPS"], yHAPS = results["yHAPS"], path = path)
        println("Tex file create in file ", path)
    else
        println("Invalid dictionary argument. It should contain the keys \"instancePath\" and \"isClientCovered\". It should also contain either \"isSiteOpened\" or (\"xHAPS\" and \"yHAPS\")")
    end 
end 

"""
Create a tex file which represents the instance and optionally a solution
"""
function texDocument(instance::Instance; clientCoverage::Vector{Float64}=Vector{Float64}([]), path::String="../results/output.tex", siteOpening::Vector{Float64}=Vector{Float64}([]), xHAPS::Vector{Float64}=Vector{Float64}([]), yHAPS::Vector{Float64}=Vector{Float64}([]))

    # Create the output folder if necessary
    mkpath(dirname(path))

    if !isfile(path)
        touch(path)
    end 

    displayASolution = length(clientCoverage) != 0

    open(path, "w") do f
        write(f, "% date = $(now()) \n")
        write(f, "\\documentclass{article} \n \\usepackage{tikz} \n \\usepackage[a4paper, total={20cm, 28cm}]{geometry}  \n \\begin{document}")
    end

    if displayASolution
        open(path, "a") do f
            write(f, "\n\\section{Solution}")
        end
    end 
    
    # Add the picture of the solution
    if length(siteOpening) != 0
        texPicture(instance, clientCoverage=clientCoverage, siteOpening = siteOpening, path = path)
    elseif length(xHAPS) != 0
        texPicture(instance, clientCoverage=clientCoverage, xHAPS = xHAPS, yHAPS = yHAPS, path = path)
    end

    open(path, "a") do f
        write(f, "\\newpage \n\n")
        if displayASolution
            write(f, "\n\\section{Instance}\n")
        end
    end 
    
    # Add the picture of the instance
    texPicture(instance, path = path)

    open(path, "a") do f
        write(f, "\\end{document} \n")
    end 
end 

"""
Display an instance or one of its solutions

Input
- instance: the instance
- (optional) clientCoverage[i] is 1 iff client i is covered (nothing if the instance is displayed rather than a solution)
- (optional) siteOpening[j] is 1 iff site j is opened (nothing if the instance is displayed rather than a solution)
- (optional) path: path in which the tex file is saved
"""
function texPicture(instance::Instance; clientCoverage::Vector{Float64}=Vector{Float64}([]), siteOpening::Vector{Float64}=Vector{Float64}([]), xHAPS::Vector{Float64}=Vector{Float64}([]),  yHAPS::Vector{Float64}=Vector{Float64}([]), path::String="../results/output.tex")

    xBounds = [minimum(instance.clientPositions[:, 1]), maximum(instance.clientPositions[:, 1])]
    yBounds = [minimum(instance.clientPositions[:, 2]), maximum(instance.clientPositions[:, 2])]
    
    displayASolution = length(clientCoverage) != 0
    open(path, "a") do f
        write(f, "\n\n \\begin{tikzpicture}")

        isCovered = Vector{Bool}(zeros(instance.n))
        
        # Display clients
        for i in 1:instance.n

            percentage = 100

            if displayASolution
                percentage = round(Int, round(clientCoverage[i], digits=3) * 100)
                if percentage > 90
                    isCovered[i] = true
                end 
            end

            x, y = sPosition(instance.clientPositions[i,1], instance.clientPositions[i,2], xBounds, yBounds)
                
            write(f, "\\node[fill=green!$(percentage)!black, minimum size=.5cm] (sol-$i) at ($(x),$(y)){{$i}}; \n") 
        end
        
        # Display sites
        if length(siteOpening) != 0
            for j in 1:instance.m

                percentage = 100

                if displayASolution
                    percentage = round(Int, round(siteOpening[j], digits=3) * 100)
                end

                x, y = sPosition(instance.sitesPositions[j,1], instance.sitesPositions[j,2], xBounds, yBounds) 

                write(f, "\\node[draw, circle, white, fill=blue!$(percentage)!lightgray, minimum size=.5cm] (alt-$j) at ($(x),$(y)){{$j}}; \n")
            end

            if displayASolution
                # Display cover links
                for i in 1:instance.n
                    for j in 1:instance.m
                        if instance.d[i,j] <= instance.rCouv
                            percentage = round(Int, round(siteOpening[j], digits=3) * 100)
                            if percentage > 0.01
                                write(f, "\\draw[blue!$(percentage)!white] (sol-$i) -- (alt-$j); \n")
                            end
                        end
                    end
                end
            end
            
            # Display communication links
            if displayASolution
                for j1 in 1:instance.m
                    for j2 in 1:instance.m
                        if j1 != j2 && instance.d_com[j1,j2] <= instance.rCom
                            percentage = round(Int, round(siteOpening[j2] * siteOpening[j1], digits=3) * 100)
                            if percentage > 0.01
                                write(f, "\\draw[blue!$(percentage)!white] (alt-$j1) -- (alt-$j2); \n")
                            end
                        end
                    end
                end
            end 
        else

            # Display the HAPS with continuous coordinates if any
            for hapsId in 1:length(xHAPS)
                x, y = sPosition(xHAPS[hapsId], yHAPS[hapsId], xBounds, yBounds)  
                write(f, "\\node[draw, circle, white, fill=blue!100!lightgray, minimum size=.5cm] (alt-$(hapsId)) at ($(x),$(y)){{$hapsId}}; \n")
            end 

            for hapsId in 1:length(xHAPS)
                x1 = xHAPS[hapsId]
                y1 = yHAPS[hapsId]
                # Display communication links
                for hapsId2 in hapsId+1:length(xHAPS)
                    x2 = xHAPS[hapsId2]
                    y2 = yHAPS[hapsId2]
                    if sqrt((x1-x2)^2+(y1-y2)^2) <= instance.rCom
                        write(f, "\\draw[blue!100!white] (alt-$hapsId) -- (alt-$hapsId2); \n") 
                    end
                end

                if displayASolution
                    # Display cover links
                    for i in 1:instance.n
                        x1 = instance.clientPositions[i, 1] 
                        y1 = instance.clientPositions[i, 2] 
                        for j in 1:length(xHAPS)
                            x2 = xHAPS[j]
                            y2 = yHAPS[j]
                            if sqrt((x1-x2)^2+(y1-y2)^2+instance.L^2) <= instance.rCouv
                                write(f, "\\draw[blue!100!white] (sol-$i) -- (alt-$j); \n")
                            end
                        end
                    end
                end
            end
        end

        scale = 19
        largestDiff = max(xBounds[2]-xBounds[1], yBounds[2]-yBounds[1]) 
        rCouv = instance.rCouv/largestDiff * scale
        rCom = instance.rCom/largestDiff * scale
        
        write(f, "\\draw[align=left]  (3.5,-1)  node{Maximal client - site communication distance:};\n")
        write(f, "\\draw[blue!100!white] (0, -1.5) -- ($(rCouv), -1.5);\n")

        write(f, "\\draw[align=left]  (3.5,-2.5)  node{Maximal site - site communication distance:};\n")
        write(f, "\\draw[blue!100!white] (0, -3) -- ($(rCom), -3);\n")

        write(f, "\\end{tikzpicture} \n\n")
    end
end


function sPosition(x, y, xBounds, yBounds)
    scale = 19
    largestDiff = max(xBounds[2]-xBounds[1], yBounds[2]-yBounds[1])
    return (x - xBounds[1])/largestDiff * scale, (y - yBounds[1])/largestDiff * scale
end 

"""
Draw simpler instance and solutions (used in articles)
"""
function simpleTexPicture(instance::Instance; clientCoverage::Vector{Float64}=Vector{Float64}([]), siteOpening::Vector{Float64}=Vector{Float64}([]), xHAPS::Vector{Float64}=Vector{Float64}(), yHAPS::Vector{Float64}=Vector{Float64}())
    
    result = "\\begin{center}\\scalebox{0.15}{ \n \\begin{tikzpicture}"
    
    # Draw the units
    for i in 1:instance.n

        # If the unit is covered
        if clientCoverage[i] > 1.0 - 1e-4	
            result *= "\\node[draw, fill=gray, minimum size=2cm] (sol-$i) at ($(instance.clientPositions[i,1])/100,$(instance.clientPositions[i, 2])/100){}; \n"	
        else	
            result *= "\\node[draw, minimum size=2cm] (sol-$i) at ($(instance.clientPositions[i, 1])/100,$(instance.clientPositions[i,2])/100){}; \n"	
        end	
    end	

    # If the HAPS are placed on a given discrete set of positions
    if length(siteOpening) >= instance.m

        # Draw each position
        for j in 1:instance.m

            # If the position has an HAPS
            if siteOpening[j] > 1.0 - 1e-4
                result *= "\\node[draw, circle, fill=gray!50!black, minimum size=2cm] (alt-$j) at ($(instance.sitesPositions[j, 1])/100,$(instance.sitesPositions[j, 2])/100){}; \n"	
                result *= "\\draw[line width=1mm, gray] ($(instance.sitesPositions[j, 1])/100,$(instance.sitesPositions[j, 2])/100) circle(15/2); \n"
            else	
                result *= "\\node[circle, fill=gray,  minimum size=0.5cm] (alt-$j) at ($(instance.sitesPositions[j, 1])/100,$(instance.sitesPositions[j, 2])/100){}; \n"	
            end	
        end

        # Draw a tree connecting all the HAPS

        # To avoid cycle, an edge is drawn if and only if it connects currently unconnected HAPS
        # siteRep is the id of the connected component of each HAPS
        siteRep = Vector{Int}(collect(1:instance.m))	

        # For each pair of positions at which there is an HAPS and which are not in the same connected component
        for j in 1:instance.m 	
            for  j1 in j+1:instance.m    	
                if siteOpening[j] > 1.0 -1e-4 && siteOpening[j1] > 1.0 - 1e-4 && instance.d_com[j1,j] <= instance.rCom && siteRep[j] != siteRep[j1]	

                    # Merge the two components and draw an edge
                    j1Rep = siteRep[j1]	
                    for j2 in 1:instance.m	
                        if siteRep[j2] == j1Rep	
                            siteRep[j2] = siteRep[j]	
                        end	
                    end 	
                    result *= "\\draw[ ultra thick] (alt-$j) -- (alt-$j1); \n"	
                end	
            end	
        end	
    end 

    # If the HAPS are not placed on a given discrete set of positions
    if length(xHAPS) >= instance.p && length(yHAPS) >= instance.p

        # Draw the discrete positions in white so that the tikzpicture has the same dimensions than the discrete solution pictures 
        for j in 1:instance.m
            result *= "\\node[circle, fill=white,  minimum size=0.5cm] (alt-$j) at ($(instance.sitesPositions[j, 1])/100,$(instance.sitesPositions[j, 2])/100){}; \n"	 
        end

        # Draw each position
        for j in 1:length(xHAPS)
            result *= "\\node[draw, circle, fill=gray!50!black, minimum size=2cm] (alt-$j) at ($(xHAPS[j])/100,$(yHAPS[j])/100){}; \n"
            result *= "\\draw[line width=1mm, gray] ($(xHAPS[j])/100,$(yHAPS[j])/100) circle(15/2); \n" 
        end

    
        # Draw a tree connecting all the HAPS

        # To avoid cycle, an edge is drawn if and only if it connects currently unconnected HAPS
        # siteRep is the id of the connected component of each HAPS
        siteRep = Vector{Int}(collect(1:length(xHAPS)))	

        # For each pair of HAPS which are not in the same connected component
        for j in 1:length(xHAPS)
            for  j1 in j+1:length(xHAPS)
                dj_j1 = sqrt((xHAPS[j]-xHAPS[j1])^2+(yHAPS[j]-yHAPS[j1])^2)
                if dj_j1 <= instance.rCom && siteRep[j] != siteRep[j1]	

                    # Merge the two components and draw an edge
                    j1Rep = siteRep[j1]	
                    for j2 in 1:length(xHAPS)
                        if siteRep[j2] == j1Rep	
                            siteRep[j2] = siteRep[j]	
                        end	
                    end 	
                    result *= "\\draw[ ultra thick] (alt-$j) -- (alt-$j1); \n"	
                end	
            end	
        end	
    end
    
    result *= "\\end{tikzpicture} \n }\\end{center} \n "	
    return result	
end	
