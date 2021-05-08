import numpy as np
import meshio
import math
import datetime
import geo_generator
from geo_generator import OPENING_TYPE, NUMBER_NODES_FRAC, VOL_ELEM_RATIO, FRAC_INITIAL_OPENING, FRAC_LENGTH, MEDIUM_HEIGHT
begin_time = datetime.datetime.now()
mesh = meshio.read("/home/keveent/PyEFVLib/meshes/msh/2D/mesh_frac.msh")

BOUNDARY_SUP = 'Sup'
BOUNDARY_INF = 'Inf'

def get_nx():
    global nx
    nx = int(NUMBER_NODES_FRAC-1)
    return nx

nx = get_nx()

def get_frac_length():
    frac_length = int(FRAC_LENGTH)
    return frac_length

def get_nodes_coord_sup(BOUNDARY_SUP):
    physical_index_sup = mesh.field_data[BOUNDARY_SUP][0]
    edges_index_sup = np.where(mesh.cell_data['gmsh:physical'][0] == physical_index_sup)[0]
    nodes_index_sup = mesh.cells[0][1][edges_index_sup]
    nodes_coord_sup = mesh.points[np.unique(nodes_index_sup)]
    return nodes_index_sup, nodes_coord_sup

def get_nodes_coord_inf(BOUNDARY_INF):
    physical_index_inf = mesh.field_data[BOUNDARY_INF][0]
    edges_index_inf = np.where(mesh.cell_data['gmsh:physical'][0] == physical_index_inf)[0]
    nodes_index_inf = mesh.cells[0][1][edges_index_inf]
    nodes_coord_inf = mesh.points[np.unique(nodes_index_inf)]
    return nodes_index_inf, nodes_coord_inf

NODES_INDEX_SUP, NODES_COORD_SUP = get_nodes_coord_sup(BOUNDARY_SUP)
NODES_INDEX_INF, NODES_COORD_INF = get_nodes_coord_inf(BOUNDARY_INF)

if nx/len(NODES_INDEX_SUP) != VOL_ELEM_RATIO or nx/len(NODES_INDEX_INF) != VOL_ELEM_RATIO:
    print('Proporção de volumes para elementos na fronteira não respeita a razão estabelecida!')

def get_x():
    # xu = np.linspace(0, FRAC_LENGTH, num=NUMBER_NODES_FRAC, endpoint=True)
    dx = geo_generator.get_dx()
    xu = np.zeros(NUMBER_NODES_FRAC)
    for i in range(1, NUMBER_NODES_FRAC):
        xu[i] =  xu[i-1] + dx[i-1]
    xp = []
    for i in range(0, nx):
        xp.append((xu[i]+xu[i+1])/2)
    xp = np.array(xp, dtype=float)
    return xp, xu

def get_width():
    xp, xu = get_x()
    if OPENING_TYPE == 'Constante':
        width = np.full((NUMBER_NODES_FRAC,), FRAC_INITIAL_OPENING) #Abertura Constante
    elif OPENING_TYPE == 'Variável': #Abertura Variável
        width = np.zeros(NUMBER_NODES_FRAC)
        width[0] = FRAC_INITIAL_OPENING/2
        width[-1] = 0.0
        slope = (width[-1]-(width[0] ))/(FRAC_LENGTH-0) #y2-y2/x2-x1
        for i in range(1, NUMBER_NODES_FRAC-1):
            width[i] = width[i-1] + ((xu[i]-xu[i-1])*slope)
        for i in range(0, len(width)):
            width[i] = width[i]*2
    width_center = []
    for i in range(0, nx):
        width_center.append((width[i]+width[i+1])/2)
    width_center = np.array(width_center, dtype=float)
    return width, width_center

def get_neighbors_sup():
    elements_neighbors_sup = []
    for i in range(0, len(NODES_INDEX_SUP)):
        for j in range(0, len(mesh.cells[-1][1])):
            if NODES_INDEX_SUP[i][0] in mesh.cells[-1][1][j] and NODES_INDEX_SUP[i][1] in mesh.cells[-1][1][j]:
                elements_neighbors_sup.append(j)
    elements_neighbors_sup = np.array(elements_neighbors_sup, dtype=int)
    return elements_neighbors_sup

def get_neighbors_inf():
    elements_neighbors_inf = []
    for i in range(0, len(NODES_INDEX_INF)):
        for j in range(0, len(mesh.cells[-1][1])):
            if NODES_INDEX_INF[i][0] in mesh.cells[-1][1][j] and NODES_INDEX_INF[i][1] in mesh.cells[-1][1][j]:
                elements_neighbors_inf.append(j)
    elements_neighbors_inf = np.array(elements_neighbors_inf, dtype=int)
    return elements_neighbors_inf

def get_nodes_neighbors_sup():
    elements_neighbors_sup = get_neighbors_sup()
    nodes_neighbors_sup = mesh.cells[-1][1][elements_neighbors_sup]
    return nodes_neighbors_sup

def get_nodes_neighbors_inf():
    elements_neighbors_inf = get_neighbors_inf()
    nodes_neighbors_inf = mesh.cells[-1][1][elements_neighbors_inf]
    return nodes_neighbors_inf

def get_centerx_neighbors():
    nodes_neighbors_sup = get_nodes_neighbors_sup()
    nodes_neighbors_inf = get_nodes_neighbors_inf()
    centerx_neighbors_sup = []
    centerx_neighbors_inf = []
    for i in range(0, len(NODES_INDEX_SUP)):
        centerx_neighbors_sup.append(np.mean(mesh.points[nodes_neighbors_sup[i]][:,0]))
        centerx_neighbors_inf.append(np.mean(mesh.points[nodes_neighbors_inf[i]][:,0]))
    centerx_neighbors_sup = np.array(centerx_neighbors_sup, dtype=float)
    centerx_neighbors_inf = np.array(centerx_neighbors_inf, dtype=float)
    return centerx_neighbors_sup, centerx_neighbors_inf

def get_centery_neighbors():
    nodes_neighbors_sup = get_nodes_neighbors_sup()
    nodes_neighbors_inf = get_nodes_neighbors_inf()
    centery_neighbors_sup = []
    centery_neighbors_inf = []
    for i in range(0, len(NODES_INDEX_SUP)):
        centery_neighbors_sup.append(np.mean(mesh.points[nodes_neighbors_sup[i]][:,1]))
        centery_neighbors_inf.append(np.mean(mesh.points[nodes_neighbors_inf[i]][:,1]))
    centery_neighbors_sup = np.array(centery_neighbors_sup, dtype=float)
    centery_neighbors_inf = np.array(centery_neighbors_inf, dtype=float)
    return centery_neighbors_sup, centery_neighbors_inf

def get_centery_faces():
    D, Dp = get_width()
    centery_sup = []
    for i in range(0, nx):
        centery_sup.append((((MEDIUM_HEIGHT+D[i])/2)+((MEDIUM_HEIGHT+D[i+1])/2))/2)
    centery_sup = np.array(centery_sup, dtype=float)
    centery_inf = []
    D = np.flip(D)
    for i in range(0, nx):
        centery_inf.append((((MEDIUM_HEIGHT-D[i])/2)+((MEDIUM_HEIGHT-D[i+1])/2))/2)
    centery_inf = np.array(centery_inf, dtype=float)
    return centery_sup, centery_inf

def get_Lns():
    centery_neighbors_sup, centery_neighbors_inf = get_centery_neighbors()
    centerx_neighbors_sup, centerx_neighbors_inf = get_centerx_neighbors()
    centery_sup, centery_inf = get_centery_faces()
    centery_neighbors_inf = np.flip(centery_neighbors_inf)
    centerx_neighbors_inf = np.flip(centerx_neighbors_inf)
    centery_inf = np.flip(centery_inf)
    xp, xu = get_x()
    Ln = []
    Ls = []
    i = 0
    k = 0
    while i < nx:
        j = 1
        while j <= VOL_ELEM_RATIO:
            Ln.append(math.hypot(xp[i] - centerx_neighbors_sup[k], centery_sup[i] - centery_neighbors_sup[k]))
            Ls.append(math.hypot(xp[i] - centerx_neighbors_inf[k], centery_inf[i] - centery_neighbors_inf[k]))
            i = i + 1
            j = j + 1
        k = k + 1
    Ln = np.array(Ln, dtype=float)
    Ls = np.array(Ls, dtype=float)
    return Ln, Ls

def get_ds():
    D, Dp = get_width()
    xp, xu = get_x()
    ds = []
    for i in range(0, nx):
        ds.append(math.hypot(xu[i+1] - xu[i], D[i+1] - D[i]))
    ds = np.array(ds, dtype=float)
    return ds

# dx = get_dx()
# D, Dp = get_width()
# xp, xu = get_x()
# Ln, Ls = get_Lns()
# ds = get_ds()

# print(datetime.datetime.now() - begin_time)