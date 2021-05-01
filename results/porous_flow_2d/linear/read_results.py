import numpy as np
import matplotlib.pyplot as plt
import csv
# import mesh_reader

RESULTS = []

with open('/home/keveent/PyEFVLib/results/porous_flow_2d/linear/Results.csv', 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        RESULTS.append(row)

X = [row[0] for row in RESULTS]
Y = [row[1] for row in RESULTS]
Z = [row[2] for row in RESULTS]
FINAL_PRESSURE = [row[-1] for row in RESULTS]

X = np.array(X[1:], dtype=float)
Y = np.array(Y[1:], dtype=float)
Z = np.array(Z[1:], dtype=float)
FINAL_PRESSURE = np.array(FINAL_PRESSURE[1:], dtype=float)
pressure = FINAL_PRESSURE
CONTORNO = plt.tricontourf(X, Y, FINAL_PRESSURE, 120)
COLOR_BAR = plt.colorbar()
COLOR_BAR.set_label('Pressão')
# plt.savefig('contorno_press_variable', dpi=1200)
plt.figure()

# plt.plot(X, Y, '.')
# plt.figure()

VALUES = np.linspace(10, 100, num=10, endpoint=True)
LEAK = np.repeat(VALUES, 2)
def get_analytical():
    L = X[55:66]
    p = []
    for i in range(0, len(L)):
        p.append((55*L[i])*(1e-3/1e-12) + 0)
    return p, L

Pf = FINAL_PRESSURE[55:66]
p, L = get_analytical()
plt.plot(L, p, label='Solução Analítica')
plt.plot(L, Pf, 'v', label='Solução Numérica')
plt.ylabel("Pressão [Pa]")
plt.xlabel("x [m]")
plt.legend()
# plt.savefig('variable_analytical', dpi=1200)
plt.figure

plt.show()


# def get_Pn(mesh, neighborsSup):
#     Pn = []
#     for i in range(0, len(neighborsSup)):
#         Pn.append((finalPressure[mesh.cells[1][1][neighborsSup[i]][0]]+finalPressure[mesh.cells[1][1][neighborsSup[i]][1]]+finalPressure[mesh.cells[1][1][neighborsSup[i]][2]]+finalPressure[mesh.cells[1][1][neighborsSup[i]][3]])/4)    
#     Pn = np.array(Pn, dtype=float)   
#     return Pn

# def get_Ps(mesh, neighborsInf):
#     Ps = []
#     for i in range(0, len(neighborsInf)):
#         Ps.append((finalPressure[mesh.cells[1][1][neighborsInf[i]][0]]+finalPressure[mesh.cells[1][1][neighborsInf[i]][1]]+finalPressure[mesh.cells[1][1][neighborsInf[i]][2]]+finalPressure[mesh.cells[1][1][neighborsInf[i]][3]])/4)
#     Ps = np.array(Ps, dtype=float)
#     return Ps