import numpy as np
import matplotlib.pyplot as plt
import csv
import mesh_reader
import json
from mesh_reader import NODES_INDEX_SUP, NODES_INDEX_INF

BOUNDARY_SUP = 'Sup'
BOUNDARY_INF = 'Inf'

mesh = mesh_reader.mesh
xp, xu = mesh_reader.get_x()
nx = mesh_reader.get_nx()
D, Dp = mesh_reader.get_width()

RESULTS_TYPE = 'Navier-Stokes'
RESULTS_TYPE = 'Poiseuille'

with open('/home/keveent/PyEFVLib/workspace/porous_flow/linear/properties.json') as f:
    properties = json.load(f)
    
DENSITY = properties['Body']['Density'] #kg/m³
VISCOSITY = properties['Body']['Viscosity'] #Pa.s

results = []

if RESULTS_TYPE == 'Navier-Stokes':
    pressure_field_file = '/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/Results_ns.csv'
elif RESULTS_TYPE == 'Poiseuille':
    pressure_field_file = '/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/Results_pois.csv'

with open(pressure_field_file, 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        results.append(row)

X = [row[0] for row in results]
Y = [row[1] for row in results]
Z = [row[2] for row in results]
FINAL_PRESSURE = [row[-1] for row in results]
# PENULT_PRESSURE = [row[-2] for row in results]

X = np.array(X[1:], dtype=float)
Y = np.array(Y[1:], dtype=float)
Z = np.array(Z[1:], dtype=float)
FINAL_PRESSURE = np.array(FINAL_PRESSURE[1:], dtype=float)
# PENULT_PRESSURE = np.array(PENULT_PRESSURE[1:], dtype=float)

np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/final_pressure_field.csv", FINAL_PRESSURE, delimiter=",")

CONTORNO = plt.tricontourf(X, Y, FINAL_PRESSURE, 120)
COLOR_BAR = plt.colorbar()
COLOR_BAR.set_label('Pressão [Pa]')
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/contorno_press', dpi=1200)
plt.figure()

NEIGHBORS_SUP = mesh_reader.get_neighbors(NODES_INDEX_SUP)
NEIGHBORS_INF = mesh_reader.get_neighbors(NODES_INDEX_INF)

def get_Pns(medium_pressure):
    for i in range(0, len(mesh.cells)):
        if mesh.cells[i][0] == 'quad':
            n = i
        elif mesh.cells[i][0] == 'triangle':
            n = i
    Pn = []
    Ps = []
    for i in range(0, len(NEIGHBORS_SUP)):
        Pn.append(np.mean(medium_pressure[mesh.cells[n][1][NEIGHBORS_SUP][i]]))
    for i in range(0, len(NEIGHBORS_INF)):
        Ps.append(np.mean(medium_pressure[mesh.cells[n][1][NEIGHBORS_INF][i]]))
    Pn = np.array(Pn, dtype=float)
    Ps = np.array(Ps, dtype=float)
    Ps = np.flip(Ps)
    return Pn, Ps

Pn, Ps = get_Pns(FINAL_PRESSURE)
np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/pn.csv", Pn, delimiter=",")
np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/ps.csv", Ps, delimiter=",")
plt.plot(xp, Pn, label='Pressão Média Superior')
plt.plot(xp, Ps, marker='.', markerfacecolor='k', markeredgecolor='k', markersize=6, linestyle='None', label='Pressão Média Inferior')
plt.xlabel("Comprimento da Fratura [m]")
# plt.xlabel("Nós dos Elementos na Vizinhança")
plt.ylabel("Pressão [Pa]")
plt.legend()
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/p_elements.pdf', dpi=1200)
plt.figure()

def get_p_tip(results):
    for i in range(0, len(mesh.cells)):
        if mesh.cells[i][0] == 'quad':
            n = i
        elif mesh.cells[i][0] == 'triangle':
            n = i
    indice = int(np.min(mesh.cells[n][1][NEIGHBORS_SUP[-1]])+1)
    p_tip = []
    for i in range(3, len(results[0])):
        p_tip.append([row[i] for row in results][indice])
    p_tip = np.array(p_tip, dtype=float)
    return p_tip

p_tip = get_p_tip(results)
plt.plot(p_tip)
plt.xlabel("Passo de Tempo")
plt.ylabel("Pressão na Ponta da Fratura [Pa]")
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/p_tip.pdf', dpi=1200)
plt.figure()

# plt.plot(p_tip)
# plt.xlabel("Tempo [s]")
# plt.ylabel("Pressão na Ponta da Fratura [Pa]")
# plt.grid()
# plt.savefig('p_tip2.pdf', dpi=1200)
# plt.figure()

# def get_face_pressure(FINAL_PRESSURE):
#     NODES_INDEX_SUP, NODES_COORD_SUP = mesh_reader.get_nodes_coord_sup(BOUNDARY_SUP)
#     NODES_INDEX_INF, NODES_COORD_INF = mesh_reader.get_nodes_coord_inf(BOUNDARY_INF)
#     face_sup_indices = np.unique(NODES_INDEX_SUP)
#     face_inf_indices = np.unique(NODES_INDEX_INF)
#     face_sup = []
#     face_inf = []
#     for i in range(0, len(face_sup_indices)):
#         face_sup.append(FINAL_PRESSURE[face_sup_indices[i]])
#     for i in range(0, len(face_inf_indices)):
#         face_inf.append(FINAL_PRESSURE[face_inf_indices[i]])
#     face_sup = np.array(face_sup, dtype=float)
#     face_inf = np.array(face_inf, dtype=float)
#     face_inf = np.flip(face_inf)
#     return face_sup, face_inf
    
    
# # face_sup, face_inf = get_face_pressure(FINAL_PRESSURE)
# # plt.ticklabel_format(axis="y", useOffset=False, style="plain")
# # plt.plot(xu, face_sup, linewidth=1.5, label='Pressão na Face Superior')
# # plt.plot(xu, face_inf, marker='.', markerfacecolor='k', markeredgecolor='k', markersize=6, linestyle='None', label='Pressão na Face Inferior')
# # # plt.xlabel("Comprimento da Fratura [m]")
# # plt.xlabel("Número de Elementos na Vizinhança")
# # plt.ylabel("Pressão [Pa]")
# # plt.legend()
# # plt.grid()
# # plt.savefig('p_faces.pdf', dpi=1200)
# # plt.figure()

plt.show()