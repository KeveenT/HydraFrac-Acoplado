from scipy.linalg import lu_factor, lu_solve
import matplotlib.pyplot as plt
import numpy as np
# import csv
import json
import mesh_reader
import leakoff
import geo_generator
from geo_generator import VOL_ELEM_RATIO

with open('/home/keveent/PyEFVLib/workspace/porous_flow/linear/properties.json') as f:
    properties = json.load(f)

with open('/home/keveent/PyEFVLib/workspace/porous_flow/linear/NumericalSettings.json') as g:
    numericalSettings = json.load(g)

#Propriedades
DENSITY = properties['Body']['Density'] #1/Pa
FLUID_COMPRESS = properties['Body']['FluidCompressibility'] #1/Pa
DENSITY = properties['Body']['Density'] #kg/m³
VISCOSITY = properties['Body']['Viscosity'] #Pa.s
PERMEABILITY = properties['Body']['Permeability'] #m²

#Geometria
FRAC_LENGTH = mesh_reader.get_frac_length() #Comprimento da Fratura
D, Dp = mesh_reader.get_width()

#Malha e Tempo
dt = numericalSettings['TimeStep']
tf = numericalSettings['FinalTime']
nt = int(tf / dt)
tol = numericalSettings['Tolerance']

#Características de Malha
nx = mesh_reader.get_nx() #Nบmero de Volumes
nu = int(nx + 1) #Nรบmero de volumes de velocidade
m = nu + nx
dxp, dxu = geo_generator.get_dxs()
Ln, Ls = mesh_reader.get_Lns()
ds = mesh_reader.get_ds()

#Montagem de Matriz de Coeficientes
def get_Auup():
    Auup = np.zeros(nu)
    Auup[0] = ((DENSITY * dxu[0] * D[0]) / (dt)) + ((2 * VISCOSITY * Dp[0]) / (dxu[0])) + 2*((6 * dxu[0] * VISCOSITY)/ (D[0]))
    if D[-1] == 0.0:
        Auup[-1] = 1
    else:
        Auup[-1] = ((DENSITY * dxu[-1] * D[-1]) / (dt)) + ((2 * VISCOSITY * Dp[-1]) / (dxu[-1])) + 2*((6 * dxu[-1] * VISCOSITY)/ (D[-1]))
    for i in range(1, len(D)-1):
        Auup[i] = ((DENSITY * dxu[i] * D[i]) / (dt)) + ((2 * VISCOSITY * Dp[i]) / (dxu[i])) + ((2 * VISCOSITY * Dp[i-1]) / (dxu[i])) + 2*((6 * dxu[i] * VISCOSITY)/ (D[i])) 
    return Auup

def get_Auuw():
    Auuw = np.zeros(nu-1)
    for i in range(0, len(Auuw)):
        Auuw[i] = 2*(VISCOSITY * Dp[i]) / (dxp[i])
    return Auuw

def get_Auue():
    Auue = np.zeros(nu-1)
    for i in range(0, len(Auue)):
        Auue[i] = 2*(VISCOSITY * Dp[i]) / (dxp[i])
    return Auue

def get_Aupe():
    Aupe = np.zeros(nx)
    for i in range(0, len(Aupe)):
        Aupe[i] = D[i]
    return Aupe

def get_Aupw():
    Aupw = np.zeros(nx)
    for i in range(0, len(Aupw)):
        Aupw[i] = D[i]
    return Aupw

def get_Aupp():
    Aupp = np.zeros(nx)
    for i in range(0, len(Aupp)):
        Aupp[i] = ((DENSITY * ds[i] * PERMEABILITY) / VISCOSITY) * ((1/Ln[i]) + (1/Ls[i]))
    return Aupp

#Montagem de Matriz de Coeficientes
def build_A(): #Monta a matriz A
    A = np.zeros([m,m])
    Auup = get_Auup()
    Auuw = get_Auuw()
    Auue = get_Auue()
    Aupw = get_Aupw()
    Aupe = get_Aupe()
    Aupp = get_Aupp()
    for i in range(0, len(Auup)):
        A[i, i] = Auup[i]
    for i in range(0, len(Auue)):
        A[i, i+1] = -Auue[i]
    for i in range(0, len(Auuw)):
        A[i+1, i] = -Auuw[i]
    for i in range(1, nx): #Adiciona o gradiente de pressão à matriz de QM
        A[i, nx+i] = -Aupw[i]
        A[i, nx+i+1] = Aupe[i]
    for i in range(0, len(Aupp)):
        A[nu+i, nu+i] = Aupp[i]
    A[0, nu] = 2 * ((D[0] + Dp[0])/2) #Condições de Contorno de pressão prescrita
    # for i in range(1, nx+nu):
    #     A[0, i] = 0
    # A[0, 0] = 1
    for i in range(nx, nx+nu):
        A[nx, i-1] = 0
    # A[nx, nx-1] = -1
    A[nx, nx] = 1
    return A

def get_A(S_star): #Retorna a matriz A com os termos advectivos
    A = build_A()
    ue = []
    uw = []
    for i in range(0, nu-1): #Velocidade da conservação da massa 
        ue.append((S_star[i]+S_star[i+1])/2)
    for i in range(1, nu):
        uw.append((S_star[i]+S_star[i-1])/2)
    for k in range(0, len(ue)): #Algoritmo para Upwind
        if ue[k] > 0 and uw[k] > 0:
            A[k, k] += DENSITY * ue[k] *Dp[k] #Adição de fluxo de massa a submatriz da conservação de quantidade de movimento
            A[k+1, k] += -DENSITY * uw[k] * Dp[k]
        else:
            A[k, k+1] += DENSITY * ue[k] * Dp[k]
            A[k+1, k+1] += -DENSITY * uw[k] * Dp[k]
    # A[nx, nx] += DENSITY * ue[-1] *Dp[-1]
    for m in range(0, nx): #Adição dos fluxos de massa a submatriz da conservação da massa
        A[nu+m, m] += -DENSITY * D[m]
    for m in range(0, nx):
        A[nu+m, m+1] += DENSITY * D[m+1]
    A[0, nu] = 2 * ((D[0] + Dp[0])/2) #Condições de Contorno de pressão prescrita
    # for i in range(1, nx+nu):
    #     A[0, i] = 0
    # A[0, 0] = 1
    for i in range(nx, nx+nu):
        A[nx, i-1] = 0
    # A[nx, nx-1] = -1
    A[nx, nx] = 1
    return A

#Construção do Termo Fonte
def get_b(S_old, Pn, Ps):
    P_entrada = 5e2
    b = np.zeros((nu+nx,1))
    # massa = 0.005*DENSITY
    for i in range(0, nu):
        b[i] = ((DENSITY * dxu[i] * D[i] * S_old[i]) / (dt))
    # b[nx] = 0
    i = 0
    k = 0
    while i < nx:
        j = 1
        while j <= VOL_ELEM_RATIO:
            b[i+nu] += ((((DENSITY * ds[i] * PERMEABILITY) / VISCOSITY) * (1 / Ln[i])) * Pn[k]) + ((((DENSITY * ds[i] * PERMEABILITY) / VISCOSITY) * (1 / Ls[i])) * Ps[k])
            i = i + 1
            j = j + 1
        k = k + 1
    b[0] +=  2 * P_entrada * ((D[0] + Dp[0])/2) #+ ((DENSITY * (dxu[0]) * D[0] * S_old[0]) / (dt))
    # b[-1] = 0
    # b[0] += massa/(DENSITY * D[0])
    return b

#Solver
def get_initial_S():
    S = np.zeros((nu+nx,1))
    return S

def solve(S_old, Pn, Ps):
    b = get_b(S_old, Pn, Ps)
    n = 0
    while n < 1:
        A = get_A(S_old)
        lu, piv = lu_factor(A)
        S_star = lu_solve((lu, piv), b)
        n = n + 1
    # if np.allclose(np.dot(A, S_star), b) == True:
    #     print('Solução Encontrada')
    # else:
    #     print('Não foi obtida uma solução')
    return S_star

def get_u_and_p(S_star, u, Pf):
    for i in range(0, nu):
        u.append(S_star[i])
    for i in range(nu, nu+nx):
        Pf.append(S_star[i])
    return u, Pf

def checkS(S_Old, S):
    Sfrac = S[-1]
    print('Diferença Solução:', max(abs((Sfrac - S_Old)/(max(Sfrac) - min(Sfrac))))[0])
    if max(abs((Sfrac - S_Old)/(max(Sfrac) - min(Sfrac)))) < tol:
        check = 'Converged'
        print('Pressão na Fratura Convergida')
        # if np.allclose(np.dot(A, S_star), b) == False:
        #     print('Vetor Solução não satisfaz as matrizes')
    else:
        check = 'Not Converged'
    return check