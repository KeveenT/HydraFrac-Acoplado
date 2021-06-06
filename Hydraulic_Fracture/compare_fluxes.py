import numpy as np
import matplotlib.pyplot as plt
import csv
import mesh_reader
import json
import geo_generator

with open('/home/keveent/PyEFVLib/workspace/porous_flow/linear/properties.json') as f:
    properties = json.load(f)
    
DENSITY = properties['Body']['Density'] #kg/m³
VISCOSITY = properties['Body']['Viscosity'] #Pa.s
PERMEABILITY = properties['Body']['Permeability'] 

D, Dp = mesh_reader.get_width()
xp, xu = mesh_reader.get_x()
nx = mesh_reader.get_nx()
dx = geo_generator.get_dx()
dxp, dxu = geo_generator.get_dxs()

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

def unique(array):
    uniq, index = np.unique(array, return_index=True)
    return uniq[index.argsort()]

def get_leakoff():
    # ds = mesh_reader.get_ds()
    leak_sup = []
    leak_inf = []
    with open(leak_sup_file, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            leak_sup.append(row)
    with open(leak_inf_file, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            leak_inf.append(row)
    leak_sup = np.array(leak_sup, dtype=float)
    leak_inf = np.array(leak_inf, dtype=float)
    # leak_sup = np.unique(leak_sup)
    # leak_inf = np.unique(leak_inf)
    leak_sup = unique(leak_sup)
    leak_inf = unique(leak_inf)
    leak_inf = np.flip(leak_inf)
    leakoff = np.zeros(len(leak_sup))
    for i in range(0, len(leakoff)):
        leakoff[i] = leak_sup[i] + leak_inf[i]
    return leakoff

def get_leakoff_ds(leakoff):
    ds = mesh_reader.get_ds()
    leak_sup = []
    leak_inf = []
    with open(leak_sup_file, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            leak_sup.append(row)
    with open(leak_inf_file, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            leak_inf.append(row)
    leak_sup = np.array(leak_sup, dtype=float)
    leak_inf = np.array(leak_inf, dtype=float)
    # leak_sup = np.unique(leak_sup)
    # leak_inf = np.unique(leak_inf)
    leak_sup = unique(leak_sup)
    leak_inf = unique(leak_inf)
    for i in range(0, len(leak_sup)):
        leak_sup[i] = ds[i] * leak_sup[i]
        leak_inf[i] = ds[i] * leak_inf[i]
    leak_inf = np.flip(leak_inf)
    leakoff = np.zeros(len(leak_sup))
    for i in range(0, len(leakoff)):
        leakoff[i] = leak_sup[i] + leak_inf[i]
    return leakoff

def get_mass_flow(u):
    mass_flow = []
    for i in range(0, len(D)):
        mass_flow.append((u[i]*D[i]))
    mass_flow = np.array(mass_flow, dtype=float)
    return mass_flow

def get_average(x):
    avg = []
    for i in range(0, len(x)-1):
        avg.append((x[i+1]+x[i])/2)
    avg = np.array(avg, dtype=float)
    return avg

def get_slope(x, y):
    slope = []
    for i in range(0, len(x)-1):
        slope.append((y[i+1]-y[i])/(x[i+1]-x[i]))
    slope = np.array(slope, dtype=float)
    # slope = abs(slope)
    return slope

def get_flow_diff(mass_flow):
    mass_flow_diff = []
    for i in range(0, len(mass_flow)-1):
        mass_flow_diff.append(mass_flow[i+1]-mass_flow[i])
    mass_flow_diff = np.array(mass_flow_diff, dtype=float)
    return mass_flow_diff

def get_permeabilities():
    frac_permeabity = np.zeros(len(D))
    for i in range(0, len(frac_permeabity)):
        frac_permeabity[i] = (D[i]**2)/(12*VISCOSITY)
    medium_permeability = np.full(len(D), 2*PERMEABILITY)
    return frac_permeabity, medium_permeability

u, Pf = get_u_and_p(RESULTS_TYPE)

leakoff = get_leakoff()
leakoff_acum = np.cumsum(leakoff)
leakoff_ds = get_leakoff_ds(leakoff)
leakoff_acum_ds = np.cumsum(leakoff_ds)

mass_flow = get_mass_flow(u)

mass_flow_acum = np.cumsum(mass_flow)
# mass_flow_slope = get_slope(xu, mass_flow)
mass_flow_slope = (np.diff(mass_flow))/(np.diff(xu))

plt.plot(xp, leakoff_acum_ds, label='Leakoff Acumulado')
plt.plot(xu, mass_flow, label='Vazão na Fratura')
plt.xlabel('Comprimento da Fratura [m]', fontsize=12)
plt.ylabel('Fluxos [m³/s]', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/flows_var.pdf', dpi=1200, bbox_inches='tight')
plt.figure()
# print(leakoff_acum_ds+mass_flow_avg)

plt.plot(xp, leakoff, label='Velocidade do Leakoff')
plt.plot(xp, mass_flow_slope, label='Inclinação da Vazão na Fratura')
plt.xlabel('Comprimento da Fratura [m]', fontsize=12)
plt.ylabel('Velocidades [m/s]', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/vels_var.pdf', dpi=1200, bbox_inches='tight')
plt.figure()
# print(leakoff+mass_flow_slope)

plt.plot(xu, mass_flow, label='Vazão na Fratura')
plt.plot(xp, leakoff_ds, label='Vazão de Leakoff')
plt.xlabel('Comprimento da Fratura [m]', fontsize=12)
plt.ylabel('Fluxos [m³/s]', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend()
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/flow_compar.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

sum_leakoff = np.sum(leakoff_ds)
# print(mass_flow[0]-sum_leakoff)
# print((mass_flow[0]-sum_leakoff)/mass_flow[0])

def mass_conservation(mass_flow, leakoff_ds):
    mass_conservation = []
    for i in range(0, len(mass_flow)-1):
        mass_conservation.append(mass_flow[i]-leakoff_ds[i]-mass_flow[i+1])
    return mass_conservation

mass_conservation = mass_conservation(mass_flow, leakoff_ds)
cm_norm = []
for i in range(0, len(mass_flow)-1):
    cm_norm.append(mass_conservation[i]/mass_flow[i])

plt.show()