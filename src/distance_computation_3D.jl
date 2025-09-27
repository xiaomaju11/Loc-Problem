using CSV
using DataFrames

function example()
    # Units
    n = 100
    map_width = 1000 # square 1000x1000
    clientPositions = rand(Float64, (n, 2)).*map_width

    # Positions
    m = 300
    sitesPositions = rand(Float64, (m, 2)).*map_width
    L = 10 # a ajuster selon la carte 3D utilisee

    # 3D map : import
    d = get_clients_sites_distances(clientPositions, sitesPositions, map_width, altitude_map, L)
end 


function get_point_altitude(x, y, map_width, altitude_map)

    vert_grid_step = map_width / (size(altitude_map)[1] - 1)
    horiz_grid_step = map_width / (size(altitude_map)[2] - 1)

    # retrouve les points de la grille qui encadrent (x,y)
    south_west_point = [floor(Int, x/horiz_grid_step), floor(Int, y/vert_grid_step)].+1
    north_west_point = [floor(Int, x/horiz_grid_step), ceil(Int, y/vert_grid_step)].+1
    south_east_point = [ceil(Int, x/horiz_grid_step), floor(Int, y/vert_grid_step)].+1
    north_east_point = [ceil(Int, x/horiz_grid_step), ceil(Int, y/vert_grid_step)].+1

    # println("Altitude : SE ", south_east_point[1], ",", south_east_point[2] ," : ")#, altitude_map[south_east_point[1], south_east_point[2]])
    # println("Altitude : NW ", north_west_point[1], ",", north_west_point[2] ," : ")#, altitude_map[north_west_point[1], north_west_point[2]])

    # on renvoie la moyenne des altitudes des points encadrants (il faudrait un barycentre)
    altitude = 0.25 * (altitude_map[south_east_point[2], south_east_point[1]]
                     + altitude_map[south_west_point[2], south_west_point[1]]
                     + altitude_map[north_east_point[2], north_east_point[1]]
                     + altitude_map[north_west_point[2], north_west_point[1]])
    return altitude
end

function get_clients_sites_distances(clientPositions, sitesPositions, map_width, altitude_map, L)
    n = size(clientPositions, 1)
    m = size(sitesPositions, 1)
    
    d = [0 for i in 1:n, j in 1:m]
    for i in 1:n, j in 1:m
        alt_client = get_point_altitude(clientPositions[i,1], clientPositions[i,2], map_width, altitude_map)
        d[i,j] = round(sqrt( (clientPositions[i, 1]-sitesPositions[j, 1])^2
                           + (clientPositions[i, 2]-sitesPositions[j, 2])^2
                           + (alt_client - L)^2 ))
    end
    return d
end

function getBrunoAltitudeMap()
    file_3D_map = "mt_bruno_elevation.csv"
    df = DataFrame(CSV.File(file_3D_map))
    altitude_map = Matrix(df)
    return altitude_map = altitude_map[:, 2:end]
end

function getBrunoAltitudeMapUpdated()
    file_3D_map = "mt_bruno_elevation_updated.csv"
    df = DataFrame(CSV.File(file_3D_map))
    altitude_map = Matrix(df)
    return altitude_map = altitude_map[:, 2:end]
end

function getBrunoAltitudeMapMax120()
    file_3D_map = "mt_bruno_elevation_max_120.csv"
    df = DataFrame(CSV.File(file_3D_map))
    altitude_map = Matrix(df)
    return altitude_map = altitude_map[:, 2:end]
end


