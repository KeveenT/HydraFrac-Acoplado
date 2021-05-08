import numpy as np
import matplotlib.pyplot as plt
import csv
import mesh_reader
import json
import geo_generator
import pandas as pd

mesh = mesh_reader.mesh
xp, xu = mesh_reader.get_x()
nx = mesh_reader.get_nx()
D, Dp = mesh_reader.get_width()
L = mesh_reader.get_frac_length()
dx = geo_generator.get_dx()

RESULTS_TYPE = 'Navier-Stokes'
# RESULTS_TYPE = 'Poiseuille'

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
    
with open('/home/keveent/PyEFVLib/workspace/porous_flow/linear/properties.json') as f:
    properties = json.load(f)
    
DENSITY = properties['Body']['Density'] #kg/m³
VISCOSITY = properties['Body']['Viscosity'] #Pa.s
PERMEABILITY = properties['Body']['Permeability']

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

def get_leakoff():
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
    leak_sup = np.unique(leak_sup)
    leak_inf = np.unique(leak_inf)
    for i in range(0, len(leak_sup)):
        leak_sup[i] = ds[i] * leak_sup[i]
        leak_inf[i] = ds[i] * leak_inf[i]
    leakoff = np.zeros(len(leak_sup))
    for i in range(0, len(leakoff)):
        leakoff[i] = leak_sup[i] + leak_inf[i]
    np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/leakoff.csv", leakoff, delimiter=",")
    return leak_sup, leak_inf

def get_mass_conservation_ns(u, leak_sup, leak_inf):
    mass_conservation_ns = []
    for i in range(0, len(D)-1):
        # mass_conservation_ns.append((u[i+1]*D[i+1])-(u[i]*D[i]-(leak_inf[i]+leak_sup[i])))
        mass_conservation_ns.append((u[i]*D[i])-(leak_inf[i]+leak_sup[i]+u[i+1]*D[i+1]))
    mass_conservation_ns = np.array(mass_conservation_ns, dtype=float)
    return mass_conservation_ns

def get_mass_flow(u):
    mass_flow = []
    for i in range(0, len(D)):
        mass_flow.append((u[i]*D[i]))
    mass_flow = np.array(mass_flow, dtype=float)
    np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/mass_flow.csv", mass_flow, delimiter=",")
    return mass_flow

def get_pressure_field():
    with open(pressure_field_file, newline='') as f:
        reader = csv.reader(f)
        row1 = next(reader)
    last_step = row1[-1]
    pressure_field = pd.read_csv(pressure_field_file)[last_step]
    return pressure_field

no_leak_iterations = get_no_leak_iterations()
plt.plot(no_leak_iterations)
plt.xlabel("Passos de Tempo")
plt.ylabel("Número de Iterações Leakoff")
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/iterations.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

u, Pf = get_u_and_p(RESULTS_TYPE)

plt.ticklabel_format(axis="y", useOffset=False, style="sci")
plt.plot(xp, Pf, linewidth=1.5)
plt.xlabel("Comprimento da Fratura [m]")
plt.ylabel("Pressão [Pa]")
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/p.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

Re = get_reynolds(u)

plt.plot(xu, Re, linewidth=1.5)
plt.xlabel("Comprimento da Fratura [m]")
plt.ylabel("Número de Reynolds Local")
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/re.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

plt.plot(xu, u, linewidth=1.5)
plt.xlabel("Comprimento da Fratura [m]")
plt.ylabel("Velocidade [m/s]")
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/v.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

leak_sup, leak_inf = get_leakoff()
plt.ticklabel_format(axis="y", useOffset=False, style="sci")
plt.plot(xp, leak_sup, label='Leakoff na Face Superior')
plt.plot(xp, leak_inf, marker='.', markerfacecolor='k', markeredgecolor='k', markersize=6, linestyle='None',  label='Leakoff na Face Inferior')
plt.xlabel("Comprimento da Fratura [m]")
plt.ylabel("Leakoff [m²/s]")
plt.grid()
plt.legend()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/leakoff.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

mass_conservation = get_mass_conservation_ns(u, leak_sup, leak_inf)
#plt.plot(xp, np.round(mass_conservation,7))
plt.plot(xp, mass_conservation)
plt.xlabel("Comprimento da Fratura")
plt.ylabel("Diferença de Vazão")  
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/conservação.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

mass_flow = get_mass_flow(u)
plt.plot(xu, mass_flow, linewidth=1.5)
plt.xlabel("Comprimento da Fratura")
plt.ylabel("Vazão")  
plt.grid()
plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/fluxo.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

# pressure_field = get_pressure_field()
# np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/final_pressure_field.csv", pressure_field, delimiter=",")

plt.show()

