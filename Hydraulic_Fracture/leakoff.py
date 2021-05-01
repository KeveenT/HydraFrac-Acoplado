import numpy as np
import matplotlib.pyplot as plt
import csv
import json
import mesh_reader
import geo_generator
from geo_generator import VOL_ELEM_RATIO

with open('/home/keveent/PyEFVLib/workspace/porous_flow/linear/properties.json') as f:
    properties = json.load(f)

with open('/home/keveent/PyEFVLib/workspace/porous_flow/linear/NumericalSettings.json') as g:
    numericalSettings = json.load(g)

with open('/home/keveent/PyEFVLib/workspace/porous_flow/linear/boundaryConditions/pressure.json') as h:
    boundaries = json.load(h)

#Propriedades
FLUID_COMPRESS = properties['Body']['FluidCompressibility'] #1/Pa
DENSITY = properties['Body']['Density'] #kg/m³
VISCOSITY = properties['Body']['Viscosity'] #Pa.s
PERMEABILITY = properties['Body']['Permeability'] #m²
INITIAL_PRESSURE = boundaries['InitialValue'] #Pa

nx = mesh_reader.get_nx() #Nบmero de Volumes
# FRAC_LENGTH = mesh_reader.get_frac_length()
# dx = geo_generator.get_dx()
# ds = mesh_reader.get_ds()
Ln, Ls = mesh_reader.get_Lns()
mesh = mesh_reader.mesh
NEIGHBORS_SUP = mesh_reader.get_neighbors_sup()
NEIGHBORS_INF = mesh_reader.get_neighbors_inf()

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
    Ps =  np.flip(Ps)
    return Pn, Ps

def get_leakoff_frac(Pf, Pn, Ps):
    Pfrac = Pf[-1]
    leakoff_sup = np.zeros(nx)
    leakoff_inf = np.zeros(nx)
    i = 0
    k = 0
    while i < nx:
        j = 1
        while j <= VOL_ELEM_RATIO:
            leakoff_sup[i] += (((PERMEABILITY / VISCOSITY) * ((Pfrac[i] - Pn[k]) / Ln[i])))
            leakoff_inf[i] += (((PERMEABILITY / VISCOSITY) * ((Pfrac[i] - Ps[k]) / Ls[i])))
            i = i + 1
            j = j + 1
        k = k + 1
    leak_sup = np.zeros(int(nx/VOL_ELEM_RATIO))
    leak_inf = np.zeros(int(nx/VOL_ELEM_RATIO))
    for i in range(0, len(leak_sup)):
        leak_sup[i] = np.mean(leakoff_sup[i*VOL_ELEM_RATIO:i*VOL_ELEM_RATIO+VOL_ELEM_RATIO])
        leak_inf[i] = np.mean(leakoff_inf[i*VOL_ELEM_RATIO:i*VOL_ELEM_RATIO+VOL_ELEM_RATIO])
    leak_sup = np.repeat(leak_sup, 2)
    leak_inf = np.repeat(leak_inf, 2)
    leak_inf = np.flip(leak_inf)
    return leak_sup, leak_inf

def leakoff_check(iterator, leak_sup, leak_inf, leak_sup_old, leak_inf_old):
    tol = 1e-3
    max = 1000
    print('Diferença Leakoff Sup:', np.max(abs((leak_sup- leak_sup_old)/(np.max(leak_sup) - np.min(leak_sup)))))
    print('Diferença Leakoff Inf:', np.max(abs((leak_inf- leak_inf_old)/(np.max(leak_inf) - np.min(leak_inf)))))
    if np.max(abs((leak_sup- leak_sup_old)/(np.max(leak_sup) - np.min(leak_sup)))) < tol and np.max(abs((leak_inf - leak_inf_old)/(np.max(leak_inf) - np.min(leak_inf)))) < tol:
        check = 'Converged'
    elif iterator == max:
        check = 'Converged'
    else:
        check = 'Not Converged'
    return check

def get_reynolds(u):
    D, Dp = mesh_reader.get_width()
    Re = np.zeros(len(D))
    for i in range(0, len(D)):
        Re[i] = (DENSITY * u[i] * 2 * D[i]) / VISCOSITY
    return Re

def get_mass_conservation(u, leak_inf, leak_sup):
    D, Dp = mesh_reader.get_width()
    leak_inf = np.unique(leak_inf)
    leak_sup = np.unique(leak_sup)
    mass_conservation = []
    for i in range(0, len(D)-1):
        mass_conservation.append((u[i+1]*D[i+1])-(u[i]*D[i]-(leak_inf[i]+leak_sup[i])))
    return mass_conservation