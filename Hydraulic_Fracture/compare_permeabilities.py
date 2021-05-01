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
# RESULTS_TYPE = 'Poiseuille'

def get_u_and_p(RESULTS_TYPE):
    if RESULTS_TYPE == 'Navier-Stokes':
        S = []
        nu = nx+1
        with open('navier_stokes.csv', 'r') as csv_file:
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
        with open('poiseuille_p.csv', 'r') as csv_file:
            csv_reader = csv.reader(csv_file)
            for row in csv_reader:
                Pf.append(row)
        Pf = Pf[-1]
        with open('poiseuille_u.csv', 'r') as csv_file:
            csv_reader = csv.reader(csv_file)
            for row in csv_reader:
                u.append(row)
        u = u[-1]
        for i in range(0, len(u)):
            u[i] = float(u[i][1:-1])
    u = np.array(u, dtype=float)
    Pf = np.array(Pf, dtype=float) 
    return u, Pf

def get_medium_permeability():
    medium_permeability = np.zeros(nx+1)
    for i in range(0, len(medium_permeability)):
        medium_permeability[i] = PERMEABILITY#/VISCOSITY
    return medium_permeability

def get_frac_permeability():
    frac_permeability = np.zeros(nx+1)
    for i in range(0, len(frac_permeability)):
        frac_permeability[i] = (D[i]**2)#/(VISCOSITY*12)
    return frac_permeability

def get_frac_gradient(RESULTS_TYPE):
    u, Pf = get_u_and_p(RESULTS_TYPE)
    frac_gradient = np.zeros(nx-1)
    for i in range(0, len(frac_gradient)-1):
        frac_gradient[i] = (-Pf[i+1]+Pf[i])/dx[i]
    return frac_gradient

def get_vol():
    frac_vol = np.zeros(len(Dp))
    for i in range(0, len(Dp)):
        frac_vol[i] = Dp[i]*dxp[i]
    return frac_vol

def get_leakoff():
    leakoff = []
    with open('leakoff.csv', 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            leakoff.append(row)
    leakoff = np.array(leakoff, dtype=float)
    return leakoff

def get_mass_flow():
    mass_flow = []
    with open('mass_flow.csv', 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            mass_flow.append(row)
    mass_flow = np.array(mass_flow, dtype=float)
    for i in range(0, len(mass_flow)):
        mass_flow[i] = mass_flow[i] * D[i]
    flow = np.zeros(int(len(mass_flow)-1))
    for i in range(0, len(flow)):
        flow[i] = (mass_flow[i]+mass_flow[i+1])/2
    return flow

medium_permeability = get_medium_permeability()
frac_permeability = get_frac_permeability()
u, Pf = get_u_and_p(RESULTS_TYPE)

plt.plot(xu, medium_permeability, label='Permeabilidade do Meio')
plt.plot(xu, frac_permeability, label='Permeabilidade da Fratura')
plt.xlabel('Comprimento da Fratura [m]')
plt.ylabel('Permeabilidade [m²]')
plt.legend()
plt.grid()
plt.figure()

xu = xu[1:-1]
frac_gradient = get_frac_gradient(RESULTS_TYPE)
plt.plot(xu, frac_gradient)
plt.xlabel('Comprimento da Fratura [m]')
plt.ylabel('Gradiente de Pressão [Pa/m]')
plt.grid()
plt.figure()

leakoff = get_leakoff()
mass_flow = get_mass_flow()
plt.plot(xp, leakoff, label='Leakoff')
plt.plot(xp, mass_flow, label='Fluxo de Massa')
plt.xlabel('Comprimento da Fratura [m]')
plt.ylabel('Fluxo [m³/s]')
plt.legend()
plt.grid()
plt.figure()

plt.show()