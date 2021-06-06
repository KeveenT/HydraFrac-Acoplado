#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 18:11:10 2021

@author: keveent
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import mesh_reader
import json
from mesh_reader import NODES_INDEX_SUP, NODES_INDEX_INF
import geo_generator

BOUNDARY_SUP = 'Sup'
BOUNDARY_INF = 'Inf'

mesh = mesh_reader.mesh
xp, xu = mesh_reader.get_x()
nx = mesh_reader.get_nx()
D, Dp = mesh_reader.get_width()
L = mesh_reader.get_frac_length()
dx = geo_generator.get_dx()
ds = mesh_reader.get_ds()
Ln, Ls = mesh_reader.get_Lns()

RESULTS_TYPE = 'Navier-Stokes'
RESULTS_TYPE = 'Poiseuille'

with open('/home/keveent/PyEFVLib/workspace/porous_flow/linear/properties.json') as f:
    properties = json.load(f)
    
DENSITY = properties['Body']['Density'] #kg/m³
VISCOSITY = properties['Body']['Viscosity'] #Pa.s
PERMEABILITY = properties['Body']['Permeability']

results = []

if RESULTS_TYPE == 'Navier-Stokes':
    pressure_field_file = '/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/Results_ns.csv'
elif RESULTS_TYPE == 'Poiseuille':
    pressure_field_file = '/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/Results_pois.csv'

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

with open(pressure_field_file, 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        results.append(row)

X = [row[0] for row in results]
Y = [row[1] for row in results]
Z = [row[2] for row in results]
FINAL_PRESSURE = [row[-1] for row in results]
FINAL_PRESSURE = np.array(FINAL_PRESSURE[1:], dtype=float)

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

def get_mass_flow(u):
    mass_flow = []
    for i in range(0, len(D)):
        mass_flow.append((u[i]*D[i]))
    mass_flow = np.array(mass_flow, dtype=float)
    return mass_flow

def get_leakoff(Pn, Ps, Pf):
    leakoff_n = []
    leakoff_s = []
    for i in range(0, len(Pf)):
        leakoff_n.append((PERMEABILITY/VISCOSITY) * ((Pf[i]-Pn[i])/Ln[i]) * ds[i]) 
        leakoff_s.append((PERMEABILITY/VISCOSITY) * ((Pf[i]-Ps[i])/Ls[i]) * ds[i])
    leakoff_n = np.array(leakoff_n, dtype=float)
    leakoff_s = np.array(leakoff_s, dtype=float)
    return leakoff_n, leakoff_s
        
Pn, Ps = get_Pns(FINAL_PRESSURE)
u, Pf = get_u_and_p(RESULTS_TYPE)
mass_flow = get_mass_flow(u)
leakoff_n, leakoff_s = get_leakoff(Pn, Ps, Pf)

def get_mass_conservation(mass_flow, leakoff_n, leakoff_s):
    mass_conservation = []
    for i in range(0, len(leakoff_n)):
        mass_conservation.append(mass_flow[i] - mass_flow[i+1] - (leakoff_n[i] + leakoff_s[i]))
    mass_conservation = np.array(mass_conservation, dtype=float)
    return mass_conservation

mass_conservation = get_mass_conservation(mass_flow, leakoff_n, leakoff_s)
plt.plot(xp, mass_conservation)
plt.xlabel("Comprimento da Fratura [m]", fontsize=12)
plt.ylabel("Diferença de Vazão", fontsize=12)  
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.grid()
# plt.savefig('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/conservação.pdf', dpi=1200, bbox_inches='tight')
plt.figure()

leakoff = leakoff_n+leakoff_s
plt.plot(xp, leakoff_n)
plt.grid()
plt.figure()

plt.plot(xp, Pn, 'o', label="Pressão na face superior")
plt.plot(xp, Pf, 'o', label="Pressão na fratura")
plt.xlabel("Comprimento da Fratura [m]", fontsize=12)
plt.ylabel("Pressão [Pa]", fontsize=12)  
plt.grid()
plt.legend()
plt.figure()