import numpy as np
import matplotlib.pyplot as plt
import csv
import json
import mesh_reader
from mesh_reader import NODES_INDEX_SUP, NODES_INDEX_INF
import geo_generator
from geo_generator import VOL_ELEM_RATIO, OPENING_TYPE

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
Ln, Ls = mesh_reader.get_Lns()
mesh = mesh_reader.mesh
NEIGHBORS_SUP = mesh_reader.get_neighbors(NODES_INDEX_SUP)
NEIGHBORS_INF = mesh_reader.get_neighbors(NODES_INDEX_INF)
if OPENING_TYPE == 'Constante':
    from mesh_reader import NODES_INDEX_TIP
    NEIGHBORS_TIP = mesh_reader.get_neighbors(NODES_INDEX_TIP)
    Ltip = mesh_reader.get_Ltip()

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

def get_Ptip(medium_pressure):
    for i in range(0, len(mesh.cells)):
        if mesh.cells[i][0] == 'quad':
            n = i
        elif mesh.cells[i][0] == 'triangle':
            n = i
    Ptip = []
    Ptip.append(np.mean(medium_pressure[mesh.cells[n][1][NEIGHBORS_TIP][0]]))
    Ptip = np.array(Ptip, dtype=float)
    return Ptip

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
    leak_sup = np.repeat(leakoff_sup, 2)
    leak_inf = np.repeat(leakoff_inf, 2)
    leak_inf = np.flip(leak_inf)
    return leak_sup, leak_inf

def get_leakoff_tip(Pf, Ptip):
    Pfrac = Pf[-1]
    leakoff_tip = np.zeros(1)
    leakoff_tip = (((PERMEABILITY / VISCOSITY) * ((Pfrac[-1] - Ptip) / Ltip)))
    leak_tip = np.repeat(leakoff_tip, 2)
    return leak_tip

def leakoff_check(iterator, leak_sup, leak_inf, leak_sup_old, leak_inf_old):
    tol = 1e-3
    max = 1000
    print('Diferença Leakoff Sup:', np.max(abs((leak_sup- leak_sup_old)/(np.max(leak_sup) - np.min(leak_sup)))))
    # print('Diferença Leakoff Inf:', np.max(abs((leak_inf- leak_inf_old)/(np.max(leak_inf) - np.min(leak_inf)))))
    if np.max(abs((leak_sup- leak_sup_old)/(np.max(leak_sup) - np.min(leak_sup)))) < tol and np.max(abs((leak_inf - leak_inf_old)/(np.max(leak_inf) - np.min(leak_inf)))) < tol:
        check = 'Converged'
    elif iterator == max:
        check = 'Converged'
    else:
        check = 'Not Converged'
    return check