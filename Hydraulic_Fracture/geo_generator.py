import numpy as np
# import matplotlib.pyplot as plt
# import meshio

#Geometry Dimensions#
FRAC_LENGTH = 10
FRAC_INITIAL_OPENING = 0.001
OPENING_TYPE = 'Variável'
TIP_OPENING = 0.0
if OPENING_TYPE == 'Constante':
    TIP_OPENING = FRAC_INITIAL_OPENING

MEDIUM_HEIGHT = 12
MEDIUM_WIDTH = 18

#Mesh Characteristics#
MESH_TYPE = 'Structured'
# MESH_TYPE = 'Unstructured'

NUMBER_NODES_FRAC = int(11)
NUMBER_VOL_FRAC = int(NUMBER_NODES_FRAC-1)

#Malha com não coincidência entre volumes e elementos
VOL_ELEM_RATIO = int(1) #Número de volumes na fratura para cada elemento em suas fronteiras
NUMBER_ELEM_NEIGHBORS = int(int(NUMBER_NODES_FRAC-1)/VOL_ELEM_RATIO)
NUMBER_NODES_NEIGHBORS = int(NUMBER_ELEM_NEIGHBORS+1)


NUMBER_NODES_SOUTH_WEST = NUMBER_NODES_NEIGHBORS#int(MEDIUM_WIDTH*(1/2))
NUMBER_NODES_NORTH_WEST = NUMBER_NODES_NEIGHBORS#int(MEDIUM_WIDTH*(1/2))

element_dx = FRAC_LENGTH/NUMBER_NODES_NORTH_WEST
no = int((MEDIUM_WIDTH-FRAC_LENGTH)/element_dx)
# no = int(1*(MEDIUM_WIDTH-FRAC_LENGTH))+1

NUMBER_NODES_SOUTH_EAST = no#int(MEDIUM_WIDTH*(1))
NUMBER_NODES_NORTH_EAST = no#int(MEDIUM_WIDTH*(1))

# no2 = int((MEDIUM_HEIGHT)/(2*element_dx))
no2 = int(4*(MEDIUM_HEIGHT))+1

NUMBER_NODES_EAST = no2#int(MEDIUM_HEIGHT*(2)+1)
NUMBER_NODES_EAST_SUP = int(NUMBER_NODES_EAST/2)
NUMBER_NODES_EAST_INF = int(NUMBER_NODES_EAST/2)
NUMBER_NODES_WEST_SUP = int((1)*NUMBER_NODES_EAST/2)
NUMBER_NODES_WEST_INF = int((1)*NUMBER_NODES_EAST/2)

NUMBER_SUPP_VERT_SUP = int(NUMBER_NODES_EAST/2)
NUMBER_SUPP_VERT_INF = int(NUMBER_NODES_EAST/2)
NUMBER_SUPP_HORIZ = NUMBER_NODES_SOUTH_EAST


def get_dx():
    common_ratio = 1
    n = NUMBER_VOL_FRAC
    if common_ratio == 1:
        dx = np.full(n, FRAC_LENGTH/n)
    else:
        summation = FRAC_LENGTH
        first_delta = (summation * (1 - common_ratio)) / (1-(common_ratio**n))
        dx = np.zeros(n)
        dx[0] = first_delta
        for i in range(1, n):
            dx[i] = dx[i-1] * common_ratio
    return dx

def get_dxs():
    dxp = get_dx()
    dxu = np.zeros(NUMBER_NODES_FRAC)
    dxu[0] = dxp[0]/2
    dxu[-1] = dxp[-1]/2
    for i in range(1, NUMBER_NODES_FRAC-1):
        dxu[i] = (dxp[i]+dxp[i-1])/2
    return dxp, dxu

def remove_repeated_rows(array, axis="row"):
    if axis == "row":
        _, indices = np.unique(array, axis=0, return_index=True)
    elif axis == "column":
        _, indices = np.unique(array, axis=1, return_index=True)
    array = array[np.sort(indices)]
    return array

def get_facets(start, nodes, loop="No"):
    if loop == "Yes":
        number_facets = int((len(nodes)))
        facets = np.zeros((number_facets, 2))
        for i in range(0, number_facets-1):
            facets[i][0] = i
            facets[i][1] = i+1
        facets[-1][0] = facets[-2][1]
        facets[-1][1] = facets[0][0]
    elif loop == "No":
        number_facets = int((len(nodes))-1)
        facets = np.zeros((number_facets, 2))
        for i in range(0, number_facets):
            facets[i][0] = i
            facets[i][1] = i+1
    facets = np.array(facets, dtype=int)
    for i in range(0, len(facets)):
        facets[i] = facets[i] + int(start)
    return facets

def add_indices(start, array):
    indices = np.zeros(len(array))
    for i in range(0, len(array)):
        indices[i] = i + start
    array = np.column_stack((indices, array))
    return array

def indices_to_string(facets):
    facets = facets.flatten() 
    facets = remove_repeated_rows(facets)
    indices = [int(i) for i in facets]
    indices = ", ".join(str(e) for e in indices)
    return indices

def get_geometry_nodes():
    global SOUTH, EAST, NORTH, WEST, FRAC
    SOUTH = np.array([[0, 0, 0], [FRAC_LENGTH, 0, 0], [MEDIUM_WIDTH, 0, 0]])
    EAST = np.array([[MEDIUM_WIDTH, 0, 0], [MEDIUM_WIDTH, (MEDIUM_HEIGHT+TIP_OPENING)/2, 0], 
                    [MEDIUM_WIDTH, (MEDIUM_HEIGHT-TIP_OPENING)/2, 0], [MEDIUM_WIDTH, MEDIUM_HEIGHT, 0]])
    EAST = remove_repeated_rows(EAST)
    NORTH = np.array([[MEDIUM_WIDTH, MEDIUM_HEIGHT, 0], [FRAC_LENGTH, MEDIUM_HEIGHT, 0], [0, MEDIUM_HEIGHT, 0]])
    WEST =  np.array([[0, MEDIUM_HEIGHT, 0], [0, (MEDIUM_HEIGHT+FRAC_INITIAL_OPENING)/2, 0], 
                     [0, (MEDIUM_HEIGHT-FRAC_INITIAL_OPENING)/2, 0], [0, 0, 0]])
    FRAC = np.array([[0, (MEDIUM_HEIGHT+FRAC_INITIAL_OPENING)/2, 0], [FRAC_LENGTH, (MEDIUM_HEIGHT+TIP_OPENING)/2, 0],
                     [FRAC_LENGTH, (MEDIUM_HEIGHT-TIP_OPENING)/2, 0], [0, (MEDIUM_HEIGHT-FRAC_INITIAL_OPENING)/2, 0]])
    FRAC = remove_repeated_rows(FRAC)
    geometry_nodes = np.concatenate((SOUTH, EAST, NORTH, WEST[0:2], FRAC, WEST[1:]))
    geometry_nodes = remove_repeated_rows(geometry_nodes)
    return geometry_nodes

def get_geometry_facets():
    geometry_facets = get_facets(1, get_geometry_nodes(), "Yes")
    return geometry_facets

def get_physicial_facets():
    south_facets = get_facets(1, SOUTH)
    east_facets = get_facets(south_facets[-1][1], EAST)
    north_facets = get_facets(east_facets[-1][1], NORTH)
    west_facets_sup = get_facets(north_facets[-1][1], WEST[0:2])
    frac_facets_sup = get_facets(west_facets_sup[-1][1], FRAC[0:2])
    if TIP_OPENING != 0.0:
        tip = get_facets(frac_facets_sup[-1][1], FRAC[1:3])
        frac_facets_inf = get_facets(tip[-1][1], FRAC[2:])
    else:
        tip = 0
        frac_facets_inf = get_facets(frac_facets_sup[-1][1], FRAC[1:])
    west_facets_inf = get_facets(frac_facets_inf[-1][1], WEST[2:])
    west_facets = np.concatenate((west_facets_sup, west_facets_inf))
    west_facets[-1][1] = 1
    return south_facets, east_facets, north_facets, west_facets, frac_facets_sup, tip, frac_facets_inf    

def get_supporting_facets():
    south_facets, east_facets, north_facets, west_facets, frac_facets_sup, tip, frac_facets_inf  = get_physicial_facets()
    vertical_suppoting_inf = np.array([[south_facets[0][1], frac_facets_inf[0][0]]])
    vertical_suppoting_sup = np.array([[frac_facets_sup[0][1], north_facets[0][1]]])
    horizontal_supporting = np.array([[frac_facets_inf[0][0], east_facets[0][1]]])
    vertical_suppoting = np.concatenate((vertical_suppoting_inf, vertical_suppoting_sup))
    return vertical_suppoting, horizontal_supporting

def get_facets_indices():
    geometry_facets = get_geometry_facets()
    south_facets, east_facets, north_facets, west_facets, frac_facets_sup, tip, frac_facets_inf = get_physicial_facets()
    south_indices = np.zeros(2)
    east_indices = np.zeros(2)
    north_indices = np.zeros(2)
    west_indices = np.zeros(2)
    for i in range(0, 2):
        south_indices[i] = int(np.where((geometry_facets == south_facets[i]).all(axis=1))[0][0]+1)
        east_indices[i] = int(np.where((geometry_facets == east_facets[i]).all(axis=1))[0][0]+1)
        north_indices[i] = int(np.where((geometry_facets == north_facets[i]).all(axis=1))[0][0]+1)
        west_indices[i] = int(np.where((geometry_facets == west_facets[i]).all(axis=1))[0][0]+1)
    frac_sup_indices = np.array([int(np.where((geometry_facets == frac_facets_sup).all(axis=1))[0][0]+1)])
    frac_inf_indices = np.array([int(np.where((geometry_facets == frac_facets_inf).all(axis=1))[0][0]+1)])
    vertical_supporting_indices = np.array([geometry_facets[-1][0]+1, geometry_facets[-1][0]+2])
    horizontal_supporting_indices = np.array([vertical_supporting_indices[-1]+1])
    return south_indices, east_indices, north_indices, west_indices, frac_sup_indices, frac_inf_indices, vertical_supporting_indices, horizontal_supporting_indices

def get_curve_loops():
    south_indices, east_indices, north_indices, west_indices, frac_sup_indices, frac_inf_indices, vertical_supporting_indices, horizontal_supporting_indices = get_facets_indices()
    south_west_quadrant = np.array([south_indices[0], vertical_supporting_indices[0],
                                    frac_inf_indices[0], west_indices[-1]])
    south_east_quadrant = np.array([south_indices[-1], east_indices[0],
                                    horizontal_supporting_indices[0], vertical_supporting_indices[0]])
    north_east_quadrant = np.array([horizontal_supporting_indices[0], east_indices[-1],
                                    north_indices[0], vertical_supporting_indices[-1]])
    north_west_quadrant = np.array([frac_sup_indices[0], vertical_supporting_indices[-1],
                                    north_indices[-1], west_indices[0]])
    return south_west_quadrant, south_east_quadrant, north_east_quadrant, north_west_quadrant

def write_mesh():
    geometry_nodes = get_geometry_nodes()
    geometry_facets = get_facets(1, geometry_nodes, "Yes")
    south_facets, east_facets, north_facets, west_facets, frac_facets_sup, tip, frac_facets_inf = get_physicial_facets()
    vertical_supporting, horizontal_supporting = get_supporting_facets()
    south_west_quadrant, south_east_quadrant, north_east_quadrant, north_west_quadrant = get_curve_loops()
    south_indices, east_indices, north_indices, west_indices, frac_sup_indices, frac_inf_indices, vertical_supporting_indices, horizontal_supporting_indices = get_facets_indices()
    
    arquivo_geo = open("mesh_frac.geo", "w")
    for i in range(0, len(geometry_nodes)):
          arquivo_geo.write(f"Point({i+1}) = {{{geometry_nodes[i][0]}, {geometry_nodes[i][1]}, {geometry_nodes[i][2]}}};\n\n")
    arquivo_geo.close()

    arquivo_geo = open("mesh_frac.geo", "a")
    length = int(len(geometry_facets)+1)
    for i in range(0, len(geometry_facets)):
            arquivo_geo.write(f"Line({i+1}) = {{{geometry_facets[i][0]}, {geometry_facets[i][1]}}};\n\n")
    arquivo_geo.write(f"Line({vertical_supporting_indices[0]}) = {{{vertical_supporting[0][0]}, {vertical_supporting[0][1]}}};\n\n")
    arquivo_geo.write(f"Line({vertical_supporting_indices[1]}) = {{{vertical_supporting[1][0]}, {vertical_supporting[1][1]}}};\n\n")
    arquivo_geo.write(f"Line({horizontal_supporting_indices[0]}) = {{{horizontal_supporting[0][0]}, {horizontal_supporting[0][1]}}};\n\n")
    
    south_east_quadrant[-1] = -south_east_quadrant[-1]
    south_east_quadrant[-2] = -south_east_quadrant[-2]
    north_east_quadrant[-1] = -north_east_quadrant[-1]
    indices1 = indices_to_string(south_west_quadrant)
    indices2 = indices_to_string(south_east_quadrant)
    indices3 = indices_to_string(north_east_quadrant)
    indices4 = indices_to_string(north_west_quadrant)
    arquivo_geo.write("Curve Loop(1) = " "{"+indices1+"}" ";\n\n")
    arquivo_geo.write("Curve Loop(2) = " "{"+indices2+"}" ";\n\n")
    arquivo_geo.write("Curve Loop(3) = " "{"+indices3+"}" ";\n\n")
    arquivo_geo.write("Curve Loop(4) = " "{"+indices4+"}" ";\n\n")
    arquivo_geo.write("Plane Surface(1) = {1};\n")
    arquivo_geo.write("Plane Surface(2) = {2};\n")
    arquivo_geo.write("Plane Surface(3) = {3};\n")
    arquivo_geo.write("Plane Surface(4) = {4};\n")

    south_facets = indices_to_string(south_indices)
    east_facets = indices_to_string(east_indices)
    north_facets = indices_to_string(north_indices)
    west_facets = indices_to_string(west_indices)
    frac_facets_sup = indices_to_string(frac_sup_indices)
    frac_facets_inf = indices_to_string(frac_inf_indices)
    arquivo_geo.write("Physical Curve("'"South"'") = " "{"+south_facets+"}" ";\n\n")
    arquivo_geo.write("Physical Curve("'"East"'") = " "{"+east_facets+"}" ";\n\n")
    arquivo_geo.write("Physical Curve("'"North"'") = " "{"+north_facets+"}" ";\n\n")
    arquivo_geo.write("Physical Curve("'"West"'") = " "{"+west_facets+"}" ";\n\n")
    arquivo_geo.write("Physical Curve("'"Sup"'") = " "{"+frac_facets_sup+"}" ";\n\n")
    arquivo_geo.write("Physical Curve("'"Inf"'") = " "{"+frac_facets_inf+"}" ";\n\n")
    arquivo_geo.write("Physical Surface("'"Body"'") =  {1, 2, 3, 4};\n\n")

    arquivo_geo.write(f"Transfinite Curve {int(south_indices[0])} = {NUMBER_NODES_SOUTH_WEST} Using Progression 1;\n\n")
    arquivo_geo.write(f"Transfinite Curve {int(south_indices[1])} = {NUMBER_NODES_SOUTH_EAST} Using Progression 1;\n\n")
    arquivo_geo.write(f"Transfinite Curve {int(east_indices[0])} = {NUMBER_NODES_EAST_INF} Using Progression 1;\n\n")
    arquivo_geo.write(f"Transfinite Curve {int(east_indices[1])} = {NUMBER_NODES_EAST_SUP} Using Progression 1;\n\n")
    arquivo_geo.write(f"Transfinite Curve {int(north_indices[0])} = {NUMBER_NODES_NORTH_EAST} Using Progression 1;\n\n")
    arquivo_geo.write(f"Transfinite Curve {int(north_indices[1])} = {NUMBER_NODES_NORTH_WEST} Using Progression 1;\n\n")
    arquivo_geo.write(f"Transfinite Curve {int(west_indices[0])} = {NUMBER_NODES_WEST_SUP} Using Progression 1;\n\n")
    arquivo_geo.write(f"Transfinite Curve {int(west_indices[1])} = {NUMBER_NODES_WEST_INF} Using Progression 1;\n\n")
    arquivo_geo.write(f"Transfinite Curve {int(frac_sup_indices[0])} = {NUMBER_NODES_NEIGHBORS} Using Progression 1;\n\n")
    arquivo_geo.write(f"Transfinite Curve {int(frac_inf_indices[0])} = {NUMBER_NODES_NEIGHBORS} Using Progression 1;\n\n")

    arquivo_geo.write(f"Transfinite Curve {int(vertical_supporting_indices[0])} = {NUMBER_SUPP_VERT_SUP} Using Progression 1;\n\n")
    arquivo_geo.write(f"Transfinite Curve {int(vertical_supporting_indices[1])} = {NUMBER_SUPP_VERT_SUP} Using Progression 1;\n\n")
    arquivo_geo.write(f"Transfinite Curve {int(horizontal_supporting_indices[0])} = {NUMBER_SUPP_HORIZ} Using Progression 1;\n\n")
    
    if MESH_TYPE == 'Structured':
        arquivo_geo.write("Transfinite Surface {1};\n")
        arquivo_geo.write("Transfinite Surface {2};\n")
        arquivo_geo.write("Transfinite Surface {3};\n")
        arquivo_geo.write("Transfinite Surface {4};\n")
        
    arquivo_geo.close()
    return

# write_mesh()
