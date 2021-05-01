import numpy as np 

FRAC_LENGTH = 10 #[m]
FRAC_INITIAL_OPENING = 0.001/5 #[m]
MEDIUM_HEIGHT = 12 #[m]
MEDIUM_WIDTH = 12 #[m]

NUMBER_NODES_FRAC = int(21)
DELTA_FRAC = FRAC_LENGTH/(NUMBER_NODES_FRAC - 1)
NUMBER_NODES_SOUTH = int(MEDIUM_WIDTH*(1))
NUMBER_NODES_NORTH = int(MEDIUM_WIDTH*(1))
NUMBER_NODES_EAST = int(MEDIUM_HEIGHT*(1))
NUMBER_NODES_EAST_SUP = int(NUMBER_NODES_EAST/2)
NUMBER_NODES_EAST_INF = int(NUMBER_NODES_EAST/2)
NUMBER_NODES_WEST_SUP = int((2)*NUMBER_NODES_EAST/1)
NUMBER_NODES_WEST_INF = int((2)*NUMBER_NODES_EAST/1)

# OPENING_TYPE = 'Constante'
OPENING_TYPE = 'Variável'

def get_frac_nodes():
    if OPENING_TYPE == 'Constante':
        frac_opening_superior = np.full((NUMBER_NODES_FRAC,), FRAC_INITIAL_OPENING) #Abertura Constante
        frac_opening_inferior = np.full((NUMBER_NODES_FRAC,), FRAC_INITIAL_OPENING)
    elif OPENING_TYPE == 'Variável': #Abertura Variável
        frac_opening_superior = np.linspace(FRAC_INITIAL_OPENING, FRAC_INITIAL_OPENING/10, 
                                   num=NUMBER_NODES_FRAC, endpoint=True)
        frac_opening_inferior = np.linspace(FRAC_INITIAL_OPENING/10, FRAC_INITIAL_OPENING, 
                                   num=NUMBER_NODES_FRAC, endpoint=True) 
    nodes_frac_superior_y = (frac_opening_superior+MEDIUM_HEIGHT)/2
    nodes_frac_inferior_y = (MEDIUM_HEIGHT-frac_opening_inferior)/2    
    nodes_frac_superior_x = np.linspace(0, FRAC_LENGTH, num=NUMBER_NODES_FRAC, endpoint=True)
    nodes_frac_inferior_x = np.linspace(FRAC_LENGTH, 0, num=NUMBER_NODES_FRAC, endpoint=True)
    return nodes_frac_inferior_y, nodes_frac_superior_y, nodes_frac_superior_x , nodes_frac_inferior_x

def get_geometry_nodes():
    nodes_south_x = np.linspace(0.0, MEDIUM_WIDTH, num=NUMBER_NODES_SOUTH, endpoint=True)
    nodes_south_y = np.full((NUMBER_NODES_SOUTH,), 0.0)

    if OPENING_TYPE == 'Constante':
        nodes_east_inferior_y = np.linspace(0.0, (MEDIUM_HEIGHT-FRAC_INITIAL_OPENING)/2, 
                                       num = NUMBER_NODES_EAST_INF, endpoint=True)
        nodes_east_superior_y = np.linspace((MEDIUM_HEIGHT+FRAC_INITIAL_OPENING)/2, MEDIUM_HEIGHT, 
                                       num = NUMBER_NODES_EAST_SUP, endpoint=True)
    elif OPENING_TYPE == 'Variável':
        nodes_east_inferior_y = np.linspace(0, MEDIUM_HEIGHT/2, 
                                       num = NUMBER_NODES_EAST_INF, endpoint=True)
        nodes_east_superior_y = np.linspace(MEDIUM_HEIGHT/2, MEDIUM_HEIGHT, 
                                       num = NUMBER_NODES_EAST_SUP, endpoint=True)
    nodes_east_x = np.full((NUMBER_NODES_EAST,), MEDIUM_WIDTH)

    nodes_north_x = np.linspace(MEDIUM_WIDTH, 0.0, num=NUMBER_NODES_SOUTH, endpoint=True)
    nodes_north_y = np.full((NUMBER_NODES_NORTH,), MEDIUM_HEIGHT)

    nodes_west_superior_y = np.linspace(MEDIUM_HEIGHT, (MEDIUM_HEIGHT+FRAC_INITIAL_OPENING)/2, 
                                   num = NUMBER_NODES_WEST_SUP, endpoint=True)
    nodes_west_superior_x = np.full((NUMBER_NODES_WEST_SUP,), 0.0)
    nodes_west_inferior_y = np.linspace((MEDIUM_HEIGHT-FRAC_INITIAL_OPENING)/2, 0.0, 
                                   num = NUMBER_NODES_WEST_INF, endpoint=True)
    nodes_west_inferior_x = np.full((NUMBER_NODES_WEST_INF,), 0.0)

    nodes_frac_inferior_y, nodes_frac_superior_y, nodes_frac_superior_x, nodes_frac_inferior_x = get_frac_nodes()

    nodes_x = np.concatenate((nodes_south_x, nodes_east_x, nodes_north_x, 
                              nodes_west_superior_x, nodes_frac_superior_x, nodes_frac_inferior_x,
                              nodes_west_inferior_x))
    nodes_y = np.concatenate((nodes_south_y, nodes_east_inferior_y, nodes_east_superior_y, 
                              nodes_north_y, nodes_west_superior_y, nodes_frac_superior_y, 
                              nodes_frac_inferior_y, nodes_west_inferior_y))
    nodes_z = np.full((len(nodes_x)), 0.0)
    geometry_nodes = np.column_stack((nodes_x, nodes_y, nodes_z))
    _, indices = np.unique(geometry_nodes, axis=0, return_index=True) 
    geometry_nodes = geometry_nodes[np.sort(indices)]
    return geometry_nodes

def get_facets(geometry_nodes):
    number_facets = int((len(geometry_nodes)))
    facets = np.zeros((number_facets, 2))
    for i in range(0, number_facets-1):
        facets[i][0] = i
        facets[i][1] = i+1
    facets[-1][0] = facets[-2][1]
    facets[-1][1] = facets[0][0]
    facets = np.array(facets, dtype=int)
    for i in range(0, len(facets)):
        facets[i] = facets[i] + 1
    return facets

def get_indices(geometry_nodes):
    indices = np.linspace(1, len(geometry_nodes), num = len(geometry_nodes))
    indices = [int(i) for i in indices]
    indices = ", ".join(str(e) for e in indices)
    return indices

def get_physical():
    south = []
    east = []
    north = []
    west = []
    superior = []
    inferior = []
    if OPENING_TYPE == 'Variável':
        for i in range(1, NUMBER_NODES_SOUTH):
            south.append(i)
        for i in range(NUMBER_NODES_SOUTH, NUMBER_NODES_SOUTH+NUMBER_NODES_EAST-2):
            east.append(i)
        for j in range(i+1, i+NUMBER_NODES_NORTH):
            north.append(j)
        for k in range(j+1, j+NUMBER_NODES_WEST_SUP):
            west.append(k)
        for l in range(k+1, k+NUMBER_NODES_FRAC):
            superior.append(l)
        ponta = 0
        for m in range(l+1, l+NUMBER_NODES_FRAC):
            inferior.append(m)
        for n in range(m+1, m+NUMBER_NODES_WEST_INF):
            west.append(n)
    elif OPENING_TYPE == 'Constante':
        for i in range(1, NUMBER_NODES_SOUTH):
            south.append(i)
        for i in range(NUMBER_NODES_SOUTH, NUMBER_NODES_SOUTH+NUMBER_NODES_EAST-1):
            east.append(i)
        for j in range(i+1, i+NUMBER_NODES_NORTH):
            north.append(j)
        for k in range(j+1, j+NUMBER_NODES_WEST_SUP):
            west.append(k)
        for l in range(k, k+NUMBER_NODES_FRAC-1):
            superior.append(l+1)
        ponta = l+2
        for m in range(l+3, l+NUMBER_NODES_FRAC+2):
            inferior.append(m)
        for n in range(m+1, m+NUMBER_NODES_WEST_INF):
            west.append(n)
    south = [int(i) for i in south]
    south = ", ".join(str(e) for e in south)
    east = [int(i) for i in east]
    east = ", ".join(str(e) for e in east)
    north = [int(i) for i in north]
    north = ", ".join(str(e) for e in north)
    west = [int(i) for i in west]
    west = ", ".join(str(e) for e in west)
    superior = [int(i) for i in superior]
    superior = ", ".join(str(e) for e in superior)
    inferior = [int(i) for i in inferior]
    inferior = ", ".join(str(e) for e in inferior)
    ponta = str(ponta)
    return south, east, north, west, superior, inferior, ponta

def write_mesh():
    geometry_nodes = get_geometry_nodes()
    facets = get_facets(geometry_nodes)
    indices = get_indices(geometry_nodes)
    south, east, north, west, superior, inferior, ponta = get_physical()
    arquivo_geo = open("mesh_frac.geo", "w")
    for i in range(0, len(geometry_nodes)):
          arquivo_geo.write(f"Point({i+1}) = {{{geometry_nodes[i][0]}, {geometry_nodes[i][1]}, {geometry_nodes[i][2]}}};\n\n")
    arquivo_geo.close()

    arquivo_geo = open("mesh_frac.geo", "a")
    for i in range(0, len(facets)):
            arquivo_geo.write(f"Line({i+1}) = {{{facets[i][0]}, {facets[i][1]}}};\n\n")
    arquivo_geo.close()

    arquivo_geo = open("mesh_frac.geo", "a")
    arquivo_geo.write("Curve Loop(1) = " "{"+indices+"}" ";\n\n")
    arquivo_geo.close()
    
    arquivo_geo = open("mesh_frac.geo", "a")
    arquivo_geo.write("Physical Curve("'"South"'") = " "{"+south+"}" ";\n\n")
    arquivo_geo.close()

    arquivo_geo = open("mesh_frac.geo", "a")
    arquivo_geo.write("Physical Curve("'"East"'") = " "{"+east+"}" ";\n\n")
    arquivo_geo.close()    

    arquivo_geo = open("mesh_frac.geo", "a")
    arquivo_geo.write("Physical Curve("'"North"'") = " "{"+north+"}" ";\n\n")
    arquivo_geo.close()
    
    arquivo_geo = open("mesh_frac.geo", "a")
    arquivo_geo.write("Physical Curve("'"West"'") = " "{"+west+"}" ";\n\n")
    arquivo_geo.close()

    arquivo_geo = open("mesh_frac.geo", "a")
    arquivo_geo.write("Physical Curve("'"Sup"'") = " "{"+superior+"}" ";\n\n")
    arquivo_geo.close()
    
    if OPENING_TYPE == 'Constante':
        arquivo_geo = open("mesh_frac.geo", "a")
        arquivo_geo.write("Physical Curve("'"Ponta"'") = " "{"+ponta+"}" ";\n\n")
        arquivo_geo.close()
    
    arquivo_geo = open("mesh_frac.geo", "a")
    arquivo_geo.write("Physical Curve("'"Inf"'") = " "{"+inferior+"}" ";\n\n")
    arquivo_geo.close()

    arquivo_geo = open("mesh_frac.geo", "a")
    arquivo_geo.write("Plane Surface(1) = {1};\n")
    arquivo_geo.close()
    
    arquivo_geo = open("mesh_frac.geo", "a")
    arquivo_geo.write("Physical Surface("'"Body"'") =  {1};\n\n")
    arquivo_geo.close()
    
    for i in range(0, len(geometry_nodes)):
        arquivo_geo = open("mesh_frac.geo", "a")
        arquivo_geo.write(f"Point{{{i+1}}} In  Surface {{1}};\n\n")
        arquivo_geo.close()
    return

write_mesh()

