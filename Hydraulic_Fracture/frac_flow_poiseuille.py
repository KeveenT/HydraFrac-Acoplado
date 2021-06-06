from scipy.linalg import lu_factor, lu_solve
import numpy as np
import matplotlib.pyplot as plt
import csv
import json
import mesh_reader
import geo_generator
import leakoff
from geo_generator import VOL_ELEM_RATIO, OPENING_TYPE

with open('/home/keveent/PyEFVLib/workspace/porous_flow/linear/properties.json') as f:
    properties = json.load(f)

with open('/home/keveent/PyEFVLib/workspace/porous_flow/linear/NumericalSettings.json') as g:
    numericalSettings = json.load(g)

#Propriedades
INJECTION = 0.005 #m3/s
FLUID_COMPRESS = properties['Body']['FluidCompressibility'] #1/Pa
DENSITY = properties['Body']['Density'] #kg/m³
VISCOSITY = properties['Body']['Viscosity'] #Pa.s
PERMEABILITY = properties['Body']['Permeability'] #m²

#Geometria
FRAC_LENGTH = mesh_reader.get_frac_length() #Comprimento da Fratura
D, Dp = mesh_reader.get_width() #Abertura

#Características de Tempo
dt = numericalSettings['TimeStep']
tf = numericalSettings['FinalTime']
tol = numericalSettings['Tolerance']
nt = int(tf / dt)

#Características de Malha
nx = mesh_reader.get_nx() #Nบmero de Volumes
dxp, dxu = geo_generator.get_dxs()
Ln, Ls = mesh_reader.get_Lns()
ds = mesh_reader.get_ds()

if OPENING_TYPE == 'Constante':
    ds_tip = mesh_reader.get_ds_tip()
    Ltip = mesh_reader.get_Ltip()

#Cálculo de Coeficientes
def get_Ae():
    Ae = np.zeros(nx)
    for i in range(0, nx):
        Ae[i] =  (((D[i+1]**3)) / (12 * VISCOSITY * dxu[i+1])) 
    return Ae

def get_Aw():
    Aw = np.zeros(nx)
    for i in range(0, nx):
        Aw[i] = (((D[i]**3)) / (12 * VISCOSITY * dxu[i]))
    return Aw

def get_Ap0():
    Ap0 = np.zeros(nx)
    for i in range(0, nx):
        Ap0[i] = ((dxp[i] * Dp[i] * FLUID_COMPRESS) / dt)
    return Ap0

def get_An():
    An = np.zeros(nx)
    for i in range(0, nx):
        An[i] = ((PERMEABILITY * ds[i]) / (VISCOSITY * Ln[i]))
    return An

def get_As():
    As = np.zeros(nx)
    for i in range(0, nx):
        As[i] = ((PERMEABILITY * ds[i]) / (VISCOSITY * Ls[i])) 
    return As

def get_Ap():
    Ap = np.zeros(nx)
    Ae = get_Ae()
    Aw = get_Aw()
    An = get_An()
    As = get_As()
    Ap0 = get_Ap0()
    for i in range(0, nx):
        Ap[i] = Ae[i] + Aw[i] + An[i] + As[i] + Ap0[i]
    return Ap

#Construção da Matriz de Coeficientes
def get_A():
    A = np.zeros((nx, nx))
    Ap = get_Ap()
    Ae = get_Ae()
    Aw = get_Aw()
    An = get_An()
    As = get_As()
    Ap0 = get_Ap0()
    for i in range(0, nx):
         A[i,i] = Ap[i]
    for j in range(1, nx):
         A[j, j-1] = -Aw[j]
    for k in range(0, nx-1):
         A[k, k+1] = -Ae[k]
    # A[0, 0] = Ae[0] + An[0] + As[0] + Ap0[0]
    A[0, 0] = Ap0[0] + Ae[0] + 2*Aw[0] + An[0] + As[0]
    A[-1, -1] = Aw[-1] + An[-1] + As[-1] + Ap0[-1]
    # A[-1, -1] = Ap0[-1] + Aw[-1] + 2*Ae[-1] + An[-1] + As[-1]
    return A

#Construção do Termo Fonte
def get_b(Pf_Old, Pn, Ps):
    global Pbw
    b = np.zeros((nx))
    Ap0 = get_Ap0()
    Aw = get_Aw()
    Ae = get_Ae()
    Pbw = 5e6
    Pbe = 1e3
    An = get_An()
    As = get_As()
    # Ps = np.flip(Ps)
    i = 0
    k = 0
    while i < nx:
        j = 1
        while j <= VOL_ELEM_RATIO:
            b[i] += Pn[k]*An[i] + Ps[k]*As[i] + Ap0[i]*Pf_Old[i]
            i = i + 1
            j = j + 1
        k = k + 1
    # b[0] +=  Pn[0]*An[0] + Ps[0]*As[0] + INJECTION + Ap0[0]*Pf_Old[0]
    b[-1] =  Pn[-1]*An[-1] + Ps[-1]*As[-1] + 0 + Ap0[-1]*Pf_Old[-1]
    b[0] = Ap0[0]*Pf_Old[0] + Pn[0]*An[0] + Ps[0]*As[0] + 2*Aw[0]*Pbw
    # b[-1] = Ap0[-1]*Pf_Old[-1] + Pn[-1]*An[-1] + Ps[-1]*As[-1] + 2*Ae[-1]*Pbe
    return b

def get_A_tip():
    A = np.zeros((nx, nx))
    Ap = get_Ap()
    Ae = get_Ae()
    Aw = get_Aw()
    An = get_An()
    As = get_As()
    Ap0 = get_Ap0()
    for i in range(0, nx):
         A[i,i] = Ap[i]
    for j in range(1, nx):
         A[j, j-1] = -Aw[j]
    for k in range(0, nx-1):
         A[k, k+1] = -Ae[k]
    # A[0, 0] = Ae[0] + An[0] + As[0] + Ap0[0]
    A[0, 0] = Ap0[0] + Ae[0] + 2*Aw[0] + An[0] + As[0]
    A[-1, -1] = Aw[-1] + An[-1] + As[-1] + Ap0[-1] + ((PERMEABILITY * ds_tip) / (VISCOSITY * Ltip))
    # A[-1, -1] = Ap0[-1] + Aw[-1] + 2*Ae[-1] + An[-1] + As[-1]
    return A

def get_b_tip(Pf_Old, Pn, Ps, Ptip):
    global Pbw
    b = np.zeros((nx))
    Ap0 = get_Ap0()
    Aw = get_Aw()
    Ae = get_Ae()
    Pbw = 1e10
    Pbe = 1e3
    An = get_An()
    As = get_As()
    # Ps = np.flip(Ps)
    i = 0
    k = 0
    while i < nx:
        j = 1
        while j <= VOL_ELEM_RATIO:
            b[i] += Pn[k]*An[i] + Ps[k]*As[i] + Ap0[i]*Pf_Old[i]
            i = i + 1
            j = j + 1
        k = k + 1
    # b[0] +=  Pn[0]*An[0] + Ps[0]*As[0] + INJECTION + Ap0[0]*Pf_Old[0]
    b[-1] =  Pn[-1]*An[-1] + Ps[-1]*As[-1] + 0 + Ap0[-1]*Pf_Old[-1] + ((PERMEABILITY * ds_tip) / (VISCOSITY * Ltip)) * Ptip
    b[0] = Ap0[0]*Pf_Old[0] + Pn[0]*An[0] + Ps[0]*As[0] + 2*Aw[0]*Pbw 
    # b[-1] = Ap0[-1]*Pf_Old[-1] + Pn[-1]*An[-1] + Ps[-1]*As[-1] + 2*Ae[-1]*Pbe
    return b

def get_initial_Pf():
    initial_Pf = np.zeros(nx)
    return initial_Pf

def solve(Pf_Old, Pn, Ps):
    A = get_A()
    b = get_b(Pf_Old, Pn, Ps)
    lu, piv = lu_factor(A)
    Pf = lu_solve((lu, piv), b)
    # if np.allclose(np.dot(A, Pf), b) == True:
    #     print('Solução Encontrada')
    # else:
    #     print('Não foi obtida uma solução')
    return Pf

def solve_tip(Pf_Old, Pn, Ps, Ptip):
    A = get_A_tip()
    b = get_b_tip(Pf_Old, Pn, Ps, Ptip)
    lu, piv = lu_factor(A)
    Pf = lu_solve((lu, piv), b)
    # if np.allclose(np.dot(A, Pf), b) == True:
    #     print('Solução Encontrada')
    # else:
    #     print('Não foi obtida uma solução')
    return Pf

def get_u(Pf):
    Pbw = 5e6
    Pbe = 1e3
    u = np.zeros((len(D),1))
    for i in range(1, len(D)-1):
        u[i] = (((Pf[i-1]-Pf[i])/dxu[i]) * ((D[i]**2) / (12 * VISCOSITY)))
    u[0] = (((Pbw-Pf[0])/(dxu[0]/2)) * ((D[0]**2) / (12 * VISCOSITY))) #CC Pressão
    # u[-1] = (((Pf[-1]-Pbe)/(dxu[-1]/2)) * ((D[-1]**2) / (12 * VISCOSITY)))
    # u[0] = INJECTION/(D[0]) #CC Vazão
    return u

def get_u_tip(Pf, Ptip):
    Pbw = 5e6
    Pbe = 1e3
    u = np.zeros((len(D),1))
    for i in range(1, len(D)-1):
        u[i] = (((Pf[i-1]-Pf[i])/dxu[i]) * ((D[i]**2) / (12 * VISCOSITY)))
    u[0] = (((Pbw-Pf[0])/(dxu[0])) * ((D[0]**2) / (12 * VISCOSITY))) #CC Pressão
    # u[-1] = (((Pf[-1]-Ptip)/((dxu[-1])+Ltip)) * ((D[-1]**2) / (12 * VISCOSITY)))
    # u[0] = INJECTION/(D[0]) #CC Vazão
    return u

def checkConvergence(Pf_Old, Pf):
    Pfrac = Pf[-1]
    print('Diferença Solução na Fratura:', max(abs((Pfrac - Pf_Old)/(max(Pfrac) - min(Pfrac)))))
    if max(abs((Pfrac - Pf_Old)/(max(Pfrac) - min(Pfrac)))) < tol:
        check = 'Converged'
        print('Pressão na Fratura Convergida')
    else:
        check = 'Not Converged'
    return check
