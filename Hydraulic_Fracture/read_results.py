import numpy as np
import matplotlib.pyplot as plt
import csv
import mesh_reader
import json
import geo_generator
from mesh_reader import NODES_INDEX_SUP, NODES_INDEX_INF

BOUNDARY_SUP = 'Sup'
BOUNDARY_INF = 'Inf'
NEIGHBORS_SUP = mesh_reader.get_neighbors(NODES_INDEX_SUP)
NEIGHBORS_INF = mesh_reader.get_neighbors(NODES_INDEX_INF)

#Características geométricas
mesh = mesh_reader.mesh
xp, xu = mesh_reader.get_x()
nx = mesh_reader.get_nx()
D, Dp = mesh_reader.get_width()
L = mesh_reader.get_frac_length()
dx = geo_generator.get_dx()
ds = mesh_reader.get_ds()
Ln, Ls = mesh_reader.get_Lns()

#Propriedades
with open('/home/keveent/PyEFVLib/workspace/porous_flow/linear/properties.json') as f:
    properties = json.load(f)

DENSITY = properties['Body']['Density'] #kg/m³
VISCOSITY = properties['Body']['Viscosity'] #Pa.s
PERMEABILITY = properties['Body']['Permeability']

#Arquivos a serem utilizados
RESULTS_TYPE = 'Navier-Stokes'
RESULTS_TYPE = 'Poiseuille'

if RESULTS_TYPE == 'Navier-Stokes':
    leak_sup_file = '/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/leakoff_superior_ns.csv'
    leak_inf_file = '/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/leakoff_inferior_ns.csv'
    iterations_file = '/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/no_leak_iterations_ns.csv'
    pressure_field_file = '/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/Results_ns.csv'
elif RESULTS_TYPE == 'Poiseuille':
    leak_sup_file = '/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/leakoff_superior_pois.csv'
    leak_inf_file = '/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/leakoff_inferior_pois.csv'
    iterations_file = '/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/no_leak_iterations_pois.csv'
    pressure_field_file = '/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/Results_pois.csv'

results = []
with open(pressure_field_file, 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        results.append(row)

X = [row[0] for row in results]
Y = [row[1] for row in results]
Z = [row[2] for row in results]
FINAL_PRESSURE = [row[-1] for row in results]

X = np.array(X[1:], dtype=float)
Y = np.array(Y[1:], dtype=float)
Z = np.array(Z[1:], dtype=float)
FINAL_PRESSURE = np.array(FINAL_PRESSURE[1:], dtype=float)

np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/final_pressure_field.csv", FINAL_PRESSURE, delimiter=",")

CONTORNO = plt.tricontourf(X, Y, FINAL_PRESSURE, 120)
COLOR_BAR = plt.colorbar()
COLOR_BAR.set_label('Pressão [Pa]', fontsize=12)
plt.xlabel("Largura [m]", fontsize=12)
plt.ylabel("Altura [m]", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
COLOR_BAR.ax.tick_params(labelsize=12) 
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/pressure_field', dpi=1200)
plt.figure()

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

def get_no_leak_iterations():
    no_leak_iterations = []
    with open(iterations_file, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            no_leak_iterations.append(row)
    no_leak_iterations = np.array(no_leak_iterations, dtype=float)
    return no_leak_iterations

def get_u_and_p(RESULTS_TYPE):
    if RESULTS_TYPE == 'Navier-Stokes':
        S = []
        nu = nx+1
        with open('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/navier_stokes.csv', 'r') as csv_file:
            csv_reader = csv.reader(csv_file)
            for row in csv_reader:
                S.append(row)
        S = S[-1]
        Sfrac = np.zeros(len(S))
        for i in range(0, len(S)):
            Sfrac[i] = float(S[i][1:-1])
        u = []
        Pf = []
        for i in range(0, nu):
            u.append(Sfrac[i])
        for i in range(nu, nu+nx):
            Pf.append(Sfrac[i])
    elif RESULTS_TYPE == 'Poiseuille':
        u = []
        Pf = []
        with open('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/poiseuille_p.csv', 'r') as csv_file:
            csv_reader = csv.reader(csv_file)
            for row in csv_reader:
                Pf.append(row)
        Pf = Pf[-1]
        with open('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/poiseuille_u.csv', 'r') as csv_file:
            csv_reader = csv.reader(csv_file)
            for row in csv_reader:
                u.append(row)
        u = u[-1]
        for i in range(0, len(u)):
            u[i] = float(u[i][1:-1])
    u = np.array(u, dtype=float)
    Pf = np.array(Pf, dtype=float) 
    return u, Pf

def get_reynolds(u):
    Re = np.zeros(len(D))
    for i in range(0, len(D)):
        Re[i] = (DENSITY * u[i] * 2 * D[i]) / VISCOSITY
    return Re

def unique(array):
    uniq, index = np.unique(array, return_index=True)
    return uniq[index.argsort()]

# def get_leakoff():
#     ds = mesh_reader.get_ds()
#     leak_sup = []
#     leak_inf = []
#     with open(leak_sup_file, 'r') as csv_file:
#         csv_reader = csv.reader(csv_file)
#         for row in csv_reader:
#             leak_sup.append(row)
#     with open(leak_inf_file, 'r') as csv_file:
#         csv_reader = csv.reader(csv_file)
#         for row in csv_reader:
#             leak_inf.append(row)
#     leak_sup = np.array(leak_sup, dtype=float)
#     leak_inf = np.array(leak_inf, dtype=float)
#     leak_sup = unique(leak_sup)
#     leak_inf = unique(leak_inf)
#     for i in range(0, len(leak_sup)):
#         leak_sup[i] = ds[i] * leak_sup[i]
#         leak_inf[i] = ds[i] * leak_inf[i]
#     leak_inf = np.flip(leak_inf)
#     leakoff = np.zeros(len(leak_sup))
#     for i in range(0, len(leakoff)):
#         leakoff[i] = leak_sup[i] + leak_inf[i]
#     np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/leakoff.csv", leakoff, delimiter=",")
#     return leak_sup, leak_inf, leakoff

def get_leakoff(Pn, Ps, Pf, area = "Yes"):
    leakoff_sup = []
    leakoff_inf = []
    if area == "Yes":
        for i in range(0, len(Pf)):
            leakoff_sup.append((PERMEABILITY/VISCOSITY) * ((Pf[i]-Pn[i])/Ln[i]) * ds[i]) 
            leakoff_inf.append((PERMEABILITY/VISCOSITY) * ((Pf[i]-Ps[i])/Ls[i]) * ds[i])
    else:
        for i in range(0, len(Pf)):
            leakoff_sup.append((PERMEABILITY/VISCOSITY) * ((Pf[i]-Pn[i])/Ln[i]) * 1) 
            leakoff_inf.append((PERMEABILITY/VISCOSITY) * ((Pf[i]-Ps[i])/Ls[i]) * 1)
    leakoff_sup = np.array(leakoff_sup, dtype=float)
    leakoff_inf = np.array(leakoff_inf, dtype=float)
    leakoff = np.zeros(len(leakoff_sup))
    for i in range(0, len(leakoff)):
        leakoff[i] = leakoff_sup[i] + leakoff_inf[i]
    return leakoff_sup, leakoff_inf, leakoff

def get_mass_flow(u):
    mass_flow = []
    for i in range(0, len(D)):
        mass_flow.append((u[i]*D[i]))
    mass_flow = np.array(mass_flow, dtype=float)
    np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/mass_flow.csv", mass_flow, delimiter=",")
    return mass_flow

def get_mass_conservation(mass_flow, leakoff_sup, leakoff_inf):
    mass_conservation = []
    for i in range(0, len(leakoff_sup)):
        mass_conservation.append(mass_flow[i] - mass_flow[i+1] - (leakoff_sup[i] + leakoff_inf[i]))
    mass_conservation = np.array(mass_conservation, dtype=float)
    return mass_conservation


no_leak_iterations = get_no_leak_iterations()
plt.plot(no_leak_iterations)
plt.xlabel("Passos de Tempo", fontsize=12)
plt.ylabel("Número de Iterações Leakoff", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/iterations.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

u, Pf = get_u_and_p(RESULTS_TYPE)
plt.ticklabel_format(axis="y", useOffset=False, style="sci")
plt.plot(xp, Pf, linewidth=1.5)
plt.xlabel("Comprimento da Fratura [m]", fontsize=12)
plt.ylabel("Pressão [Pa]", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/pressure.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

plt.plot(xu[:-1], u[:-1], linewidth=1.5)
plt.xlabel("Comprimento da Fratura [m]", fontsize=12)
plt.ylabel("Velocidade [m/s]", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/velocity.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

Re = get_reynolds(u)
plt.plot(xu, Re, linewidth=1.5)
plt.xlabel("Comprimento da Fratura [m]", fontsize=12)
plt.ylabel("Número de Reynolds Local", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/reynolds.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

mass_flow = get_mass_flow(u)
plt.plot(xu, mass_flow, linewidth=1.5)
plt.xlabel("Comprimento da Fratura [m]", fontsize=12)
plt.ylabel("Vazão Volumétrica [m³/s]", fontsize=12) 
plt.xticks(fontsize=12)
plt.yticks(fontsize=12) 
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/mass_flow.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

Pn, Ps = get_Pns(FINAL_PRESSURE)
np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/pn.csv", Pn, delimiter=",")
np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/ps.csv", Ps, delimiter=",")
plt.plot(xp, Pn, label='Pressão Média Superior')
plt.plot(xp, Ps, marker='.', markerfacecolor='k', markeredgecolor='k', markersize=6, linestyle='None', label='Pressão Média Inferior')
plt.xlabel("Comprimento da Fratura [m]", fontsize=12)
plt.ylabel("Pressão [Pa]", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/p_elements.pdf', dpi=1200)
plt.figure()

leak_sup, leak_inf, leakoff = get_leakoff(Pn, Ps, Pf)
plt.ticklabel_format(axis="y", useOffset=False, style="sci")
plt.plot(xp, leak_sup, label='Leakoff na Face Superior')
plt.plot(xp, leak_inf, marker='.', markerfacecolor='k', markeredgecolor='k', markersize=6, linestyle='None',  label='Leakoff na Face Inferior')
plt.xlabel("Comprimento da Fratura [m]", fontsize=12)
plt.ylabel("Leakoff [m³/s]", fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
plt.legend()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/leakoff.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

mass_conservation = get_mass_conservation(mass_flow, leak_sup, leak_inf)
#plt.plot(xp, np.round(mass_conservation,12))
plt.plot(xp, mass_conservation)
plt.xlabel("Comprimento da Fratura", fontsize=12)
plt.ylabel("Diferença de Vazão", fontsize=12) 
plt.xticks(fontsize=12)
plt.yticks(fontsize=12) 
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/conservação.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

plt.plot(xp, Pn, label='Pressão Média Superior')
plt.plot(xp, Pf, label="Pressão na Fratura")
plt.plot(xp, Ps, marker='.', markerfacecolor='k', markeredgecolor='k', markersize=6, linestyle='None', label='Pressão Média Inferior')
plt.xlabel("Comprimento da Fratura [m]", fontsize=12)
plt.ylabel("Pressão [Pa]", fontsize=12)  
plt.grid()
plt.legend()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/press_frac_medium.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

leakoff_acum = np.cumsum(leakoff)
mass_flow_acum = np.cumsum(mass_flow)

plt.plot(xp, leakoff_acum, label='Leakoff Acumulado')
plt.plot(xu, mass_flow, label='Vazão na Fratura')
plt.xlabel('Comprimento da Fratura [m]', fontsize=12)
plt.ylabel('Fluxos [m³/s]', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/flows_var.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

leak_sup_v, leak_inf_v, leakoff_v = get_leakoff(Pn, Ps, Pf, area = "No")
leakoff_acum_v = np.cumsum(leakoff_v)
mass_flow_slope = (np.diff(mass_flow))/(np.diff(xu))

plt.plot(xp, leakoff_v, label='Velocidade do Leakoff')
plt.plot(xp, mass_flow_slope, label='Inclinação da Vazão na Fratura')
plt.xlabel('Comprimento da Fratura [m]', fontsize=12)
plt.ylabel('Velocidades [m/s]', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/vels_var.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

plt.plot(xu, mass_flow, label='Vazão na Fratura')
plt.plot(xp, leakoff, label='Vazão de Leakoff')
plt.xlabel('Comprimento da Fratura [m]', fontsize=12)
plt.ylabel('Fluxos [m³/s]', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/flow_compar.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

plt.show()