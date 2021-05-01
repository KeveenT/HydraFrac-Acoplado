import numpy as np
import matplotlib.pyplot as plt

#Geometria e Malha da Fratura
initialFracOpening = 1 #[m]
fracLength = 10 #[m]
dy = initialFracOpening 
nxFrac = 10 #Número de volumes
dx = fracLength/nxFrac

#Geometria e Malha do Meio
mediumWidth = 15 #[m]
mediumHeight = 15 #[m]
nxMedium = mediumWidth/dx #Número de elementos
nyMedium = mediumHeight/dy #Número de elementos
if nxMedium % 1 != 0:
    print("Com essa abertura o número de elementos nas fronteiras North e South não são inteiros!")
else:
    nxMedium = int(nxMedium)
if nxMedium % 1 != 0:
    print("Com essa abertura o número de elementos nas fronteiras East e West não são inteiros!")
else:
    nyMedium = int(nyMedium)

geomPoints = np.array([[0.0, 0.0, 0.0],
                       [mediumWidth, 0.0, 0.0],
                       [mediumWidth, mediumHeight, 0.0],
                       [0.0, mediumHeight, 0.0]]) 

def get_nodesCoordinates():
  
    interpolX = np.zeros(nxMedium+1)
    for i in range(0, len(interpolX)):
        interpolX[i] = geomPoints[0][0] + dx*i
        
    interpolY = np.zeros(nyMedium+1)
    for i in range(0, len(interpolY)):
        interpolY[i] = geomPoints[1][1] + dy*i
    
    x = np.zeros(len(interpolX)*len(interpolY))
    y = np.zeros(len(interpolX)*len(interpolY))
    z = np.zeros(len(interpolX)*len(interpolY))

    i = len(interpolX)
    while i != len(x)+len(interpolX):
        j = 0
        while j != len(interpolX):
            x[i-len(interpolX)] = interpolX[j]
            j = j+1
            i = i + 1

    for i in range(0, int(len(y))):
        y[i] = interpolY[int(round(i/(nxMedium+1),5))]
    
    nodesCoordinates = np.zeros((len(x), 3))
    for i in range(0, len(nodesCoordinates)):
        nodesCoordinates[i][0] = x[i]
        nodesCoordinates[i][1] = y[i]
        nodesCoordinates[i][2] = z[i]
    
    nodesNumber = len(nodesCoordinates)
    
    return nodesCoordinates, nodesNumber

def get_facetsConectivity():
    facetsConectivity = np.zeros((int(2*(nxMedium + nyMedium + nxFrac)+1), 2)) #Esquerda para a direita, baixo para cima
    nodesCoordinates, nodesNumber = get_nodesCoordinates()
    
    #Fronteira South
    indicesSouth = np.where([row[1] for row in nodesCoordinates] == geomPoints[0][1]) 
    for i in range(0, nxMedium):
        facetsConectivity[i][0] = indicesSouth[0][i]
        facetsConectivity[i][1] = indicesSouth[0][i+1]

    #Fronteira East
    indicesEast = np.where([row[0] for row in nodesCoordinates] == geomPoints[1][0]) 
    for i in range(0, nyMedium):
        facetsConectivity[i+nxMedium][0] = indicesEast[0][i]
        facetsConectivity[i+nxMedium][1] = indicesEast[0][i+1]

    #Fronteira North 
    indicesNorth = np.where([row[1] for row in nodesCoordinates] == geomPoints[2][1]) 
    for i in range(0, nxMedium): 
        facetsConectivity[i+nxMedium+nyMedium][0] = indicesNorth[0][i]
        facetsConectivity[i+nxMedium+nyMedium][1] = indicesNorth[0][i+1]
    
    #Fronteira West
    indicesWest = np.where([row[0] for row in nodesCoordinates] == geomPoints[3][0]) 
    for i in range(0, nyMedium): 
        facetsConectivity[i+nxMedium+nyMedium+nxMedium][0] = indicesWest[0][i]
        facetsConectivity[i+nxMedium+nyMedium+nxMedium][1] = indicesWest[0][i+1]
    
    facetsConectivity = np.delete(facetsConectivity, nxMedium+nyMedium+nxMedium+int((len(indicesWest[0])-1)/2), axis=0) #Apaga a conecetividade referente a entrada da fratura

    #Fronteira Superior Fratura
    for i in range(0, nxFrac): 
        facetsConectivity[i+nxMedium+nyMedium+nxMedium+nyMedium-1][0] = facetsConectivity[nxMedium+nyMedium+nxMedium+int((len(indicesWest[0])-1)/2)][0] + i 
        facetsConectivity[i+nxMedium+nyMedium+nxMedium+nyMedium-1][1] = facetsConectivity[nxMedium+nyMedium+nxMedium+int((len(indicesWest[0])-1)/2)][0] + i + 1
    
    # Fronteira Superior Fratura
    for i in range(0, nxFrac): 
        facetsConectivity[i+nxMedium+nyMedium+nxMedium+nyMedium+nxFrac-1][0] = facetsConectivity[nxMedium+nyMedium+nxMedium+int((len(indicesWest[0])-1)/2)-1][1] + i 
        facetsConectivity[i+nxMedium+nyMedium+nxMedium+nyMedium+nxFrac-1][1] = facetsConectivity[nxMedium+nyMedium+nxMedium+int((len(indicesWest[0])-1)/2)-1][1] + i + 1

    #Ponta
    facetsConectivity[-1][0] = facetsConectivity[-2][1]
    facetsConectivity[-1][1] = facetsConectivity[-(2+nxFrac)][1]

    for i in range(0, len(facetsConectivity)):
        facetsConectivity[i][0] += 1
        facetsConectivity[i][1] += 1
    
    facetsConectivity = np.array(facetsConectivity, dtype=int)
    
    facetsNumber = len(facetsConectivity)
    return facetsConectivity, facetsNumber

def get_elementsConectivity():
    nodesCoordinates, nodesNumber = get_nodesCoordinates()
    indicesNodes = np.zeros((len(nodesCoordinates), 1))
    for i in range(0, len(indicesNodes)):
        indicesNodes[i][0] = i 

    nodesCoordinates = np.hstack((indicesNodes, nodesCoordinates))
    nodesCoordinates = nodesCoordinates[nodesCoordinates[:,1].argsort()]
    nodesCoordinates = nodesCoordinates[nodesCoordinates[:,2].argsort(kind='mergesort')]
    
    elementsConectivity = np.zeros((nxMedium*nyMedium, 4))
    
    j = 0
    for i in range(0, nxMedium*nyMedium):
        if i == 0:
            j = i
        elif i % nxMedium == 0:
            j = j + 2
        else:
            j = j + 1
        elementsConectivity[i][0] = nodesCoordinates[j][0]
        elementsConectivity[i][1] = nodesCoordinates[j+1][0]
        elementsConectivity[i][2] = nodesCoordinates[nxMedium+j+2][0]
        elementsConectivity[i][3] = nodesCoordinates[nxMedium+j+1][0]

    facetsConectivity, facetsNumber = get_facetsConectivity()
    fracConectivity = np.zeros((nxFrac*2, len(facetsConectivity[0])))
    for i in range(0, len(fracConectivity)):
        fracConectivity[i] = facetsConectivity[len(facetsConectivity)-(nxFrac*2)-1+i]
    
    elementsConectivityFrac = np.zeros((nxFrac, 4))
    for i in range(0, len(elementsConectivityFrac)):
        elementsConectivityFrac[i][0] = fracConectivity[nxFrac+i][0] 
        elementsConectivityFrac[i][1] = fracConectivity[nxFrac+i][1]
        elementsConectivityFrac[i][2] = fracConectivity[i][1]
        elementsConectivityFrac[i][3] = fracConectivity[i][0]
    
    fracElementsIndices = []
    for i in range(0, len(elementsConectivityFrac)):
        fracElementsIndices.append(np.where((elementsConectivity == elementsConectivityFrac[i])))
    
    for i in range(0, len(fracElementsIndices)):
        fracElementsIndices[i] = fracElementsIndices[i][0][0]
        fracElementsIndices[i] -= 1
    elementsConectivity = np.delete(elementsConectivity, fracElementsIndices, axis=0)

    for i in range(0, len(elementsConectivity)):
        elementsConectivity[i][0] += 1
        elementsConectivity[i][1] += 1
        elementsConectivity[i][2] += 1
        elementsConectivity[i][3] += 1        
    
    elementsNumber = len(elementsConectivity)
    
    return elementsConectivity, elementsNumber

def get_boundaryConditions():
    boundaryConditions = ['"South"', '"East"', '"North"', '"West"', '"Sup"', '"Inf"', '"Ponta"', '"Body"']
    boundaryConditionsNumber = len(boundaryConditions)
    return boundaryConditions, boundaryConditionsNumber

def write_msh():
    nodesCoordinates, nodesNumber = get_nodesCoordinates()
    # for i in range(96, 112):
    #     nodesCoordinates[i][1] = 5.75
        
    # for i in range(144, 160):
    #     nodesCoordinates[i][1] = 9.25
    facetsConectivity, facetsNumber = get_facetsConectivity()
    elementsConectivity, elementsNumber = get_elementsConectivity()
    boundaryConditions, boundaryConditionsNumber = get_boundaryConditions()
    elementsNumberTotal = facetsNumber + elementsNumber
    print('Número de Elementos:', elementsNumberTotal)
    mshFile = open("mesh_frac.msh", "w")
    mshFile.write("$MeshFormat \n2.2 0 8 \n$EndMeshFormat \n")
    mshFile.close()

    mshFile = open("mesh_frac.msh", "a")
    mshFile.write(f"$PhysicalNames \n{boundaryConditionsNumber} \n")
    for i in range(0, boundaryConditionsNumber-1):
        mshFile.write(f"1 {i+1} {boundaryConditions[i]} \n") #physical-dimension physical-number
    mshFile.write(f"2 {boundaryConditionsNumber} {boundaryConditions[-1]} \n")
    mshFile.write("$EndPhysicalNames \n")
    mshFile.close()

    mshFile = open("mesh_frac.msh", "a")
    mshFile.write(f"$Nodes \n{nodesNumber} \n")
    for i in range(0, nodesNumber):
        mshFile.write(f"{i+1} {nodesCoordinates[i][0]} {nodesCoordinates[i][1]} {nodesCoordinates[i][2]}\n")
    mshFile.write("$EndNodes \n")
    mshFile.close()

    mshFile = open("mesh_frac.msh", "a")
    mshFile.write(f"$Elements \n{elementsNumberTotal} \n")
    for i in range(0, nxMedium):
        mshFile.write(f"{i+1} 1 2 1 1 {int(facetsConectivity[i][0])} {int(facetsConectivity[i][1])}\n")    
    for i in range(nxMedium, nxMedium+nyMedium):
        mshFile.write(f"{i+1} 1 2 2 2 {int(facetsConectivity[i][0])} {int(facetsConectivity[i][1])}\n")   
    for i in range(nxMedium+nyMedium, 2*nxMedium+nyMedium):
        mshFile.write(f"{i+1} 1 2 3 3 {int(facetsConectivity[i][0])} {int(facetsConectivity[i][1])}\n")   
    for i in range(2*nxMedium+nyMedium, 2*nxMedium+nyMedium*2-1):
        mshFile.write(f"{i+1} 1 2 4 4 {int(facetsConectivity[i][0])} {int(facetsConectivity[i][1])}\n")
    for i in range(2*nxMedium+nyMedium*2-1, 2*nxMedium+nyMedium*2-1+nxFrac):
        mshFile.write(f"{i+1} 1 2 5 5 {int(facetsConectivity[i][0])} {int(facetsConectivity[i][1])}\n") 
    for i in range(2*nxMedium+nyMedium*2-1+nxFrac, 2*nxMedium+nyMedium*2-1+(2*nxFrac)):
        mshFile.write(f"{i+1} 1 2 6 6 {int(facetsConectivity[i][0])} {int(facetsConectivity[i][1])}\n")
    mshFile.write(f"{2*nxMedium+nyMedium*2-1+(2*nxFrac)+1} 1 2 7 7 {int(facetsConectivity[-1][0])} {int(facetsConectivity[-1][1])}\n") 

    p = boundaryConditions.index('"Body"') + 1
    e = 1
    for i in range(0, elementsNumber):
        mshFile.write(f"{facetsNumber+i+1} 3 2 {p} {e} {int(elementsConectivity[i][0])} {int(elementsConectivity[i][1])} {int(elementsConectivity[i][2])} {int(elementsConectivity[i][3])}\n")
    mshFile.write("$EndElements \n")
    mshFile.close()
    return

write_msh()



