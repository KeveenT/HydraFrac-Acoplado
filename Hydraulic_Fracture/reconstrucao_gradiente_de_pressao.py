import sys,os
pyEFVLibPath = os.path.join(os.path.dirname(__file__), os.path.pardir, os.path.pardir)
pyEFVLibPath = "/home/keveent/PyEFVLib/"
workspacePath = os.path.join(os.path.dirname(__file__), os.path.pardir)
sys.path += [pyEFVLibPath, workspacePath]

import PyEFVLib
import numpy as np
import mesh_reader
import json
import csv

def reconstruct2D(grid, scalarField):
    gradientField = np.zeros((len(scalarField), 2))

    for element in grid.elements:
        elementFieldVector = np.array( [scalarField[vertex.handle] for vertex in element.vertices] )
        for innerFaceIndex, innerFace in enumerate(element.innerFaces):
            backwardVertex = element.vertices[ element.shape.innerFaceNeighborVertices[innerFaceIndex][0] ]
            forwardVertex = element.vertices[ element.shape.innerFaceNeighborVertices[innerFaceIndex][1] ]

            shapeFunctionValues = element.shape.innerFaceShapeFunctionValues[innerFaceIndex]
            pressureAtIP = np.dot(elementFieldVector, shapeFunctionValues)

            area = innerFace.area.getCoordinates()[:-1]

            gradientField[backwardVertex.handle] += pressureAtIP * area    
            gradientField[forwardVertex.handle] -= pressureAtIP * area    

    for facet in grid.facets:
        elementFieldVector = np.array( [scalarField[vertex.handle] for vertex in facet.element.vertices] )
        for outerFace in facet.outerFaces:
            shapeFunctionValues = facet.element.shape.outerFaceShapeFunctionValues[facet.elementLocalIndex][outerFace.local]

            pressureAtIP = np.dot(elementFieldVector, shapeFunctionValues)
            area = outerFace.area.getCoordinates()[:-1]

            gradientField[outerFace.vertex.handle] += pressureAtIP * area    

    gradientField = np.array([ grad/vertex.volume for vertex, grad in zip( grid.vertices, gradientField ) ])

    return gradientField

if __name__ == "__main__":
    grid = PyEFVLib.Grid( PyEFVLib.MSHReader( os.path.join(pyEFVLibPath, "meshes", "msh", "2D", "mesh_frac.msh") ).getData() )

    # Importar o campo de pressões de um arquivo csv
    import pandas as pd
    fileName = "Results.csv"
    
    with open(fileName, newline='') as f:
        reader = csv.reader(f)
        row1 = next(reader)
    last_step = row1[-1]
    
    field = pd.read_csv(fileName)[last_step]
#     field = np.array([3*(vertex.x**2)+4*(vertex.y**2) for vertex in grid.vertices])

    gradientField = reconstruct2D(grid,field)
    import matplotlib.pyplot as plt
    X,Y = zip(*[vertex.getCoordinates()[:-1] for vertex in grid.vertices])
    U,V = zip(*gradientField)
    
    plt.quiver(X,Y,U,V)
    plt.savefig('gradientes.pdf', dpi=1200)
    plt.show()
    
mesh = mesh_reader.mesh

with open('/home/keveent/PyEFVLib/workspace/porous_flow/linear/properties.json') as f:
    properties = json.load(f)
    
DENSITY = properties['Body']['Density'] #kg/m³
VISCOSITY = properties['Body']['Viscosity'] #Pa.s
PERMEABILITY = properties['Body']['Permeability']

def get_face_gradients(BOUNDARY):
    NODES_INDEX, NODES_COORD = mesh_reader.get_nodes_coord_sup(BOUNDARY)
    face_indices = np.unique(NODES_INDEX)
    face_gradients = []
    for i in range(0, len(face_indices)):
        face_gradients.append(gradientField[face_indices[i]])
    face_gradients = np.array(face_gradients, dtype=float)
    return face_gradients

def get_magnitude(gradientField):
    magnitude = []
    for i in range(0, len(gradientField)):
        magnitude.append(np.sqrt((gradientField[i][0]**2)+(gradientField[i][1]**2)))
    magnitude = np.array(magnitude, dtype=float)
    magnitude_average = []
    for i in range(0, len(magnitude)-1):
        magnitude_average.append(((magnitude[i]+magnitude[i+1])/2))
    magnitude_average = np.array(magnitude_average, dtype=float)
    return magnitude_average

def get_velocity(magnitude):
    velocity = []
    for i in range(0, len(magnitude)):
        velocity.append(magnitude[i]*(PERMEABILITY/VISCOSITY))
    velocity = np.array(velocity, dtype=float)
    return velocity

def get_area(BOUNDARY):
    boundary_name = [boundary for boundary in grid.boundaries if boundary.name == BOUNDARY][0]
    i = 0
    half_area = []
    for facet in boundary_name.facets:
        for outerFace in facet.outerFaces:
            half_area.append(np.linalg.norm(outerFace.area.getCoordinates()))
            i = i + 1
    half_area = np.array(half_area, dtype=float)
    area = np.zeros(int(len(half_area)/2))
    for i in range(0, len(area)):
        area[i] = half_area[i]+half_area[i+1]
    return area

def get_flow(velocity, area):
    flow = []
    for i in range(0, len(area)):
        flow.append(velocity[i]*area[i])
    flow = np.array(flow, dtype=float)
    return flow

gradients_sup = get_face_gradients('Sup')
gradients_inf = get_face_gradients('Inf')
gradients_south = get_face_gradients('South')
gradients_east = get_face_gradients('East')
gradients_north = get_face_gradients('North')
gradients_west = get_face_gradients('West')

magnitude_sup = get_magnitude(gradients_sup)
magnitude_inf = get_magnitude(gradients_inf)
magnitude_inf = np.flip(magnitude_inf)
magnitude_south = get_magnitude(gradients_south)
magnitude_east = get_magnitude(gradients_east)
magnitude_north = get_magnitude(gradients_north)
magnitude_west = get_magnitude(gradients_west)

# np.savetxt("magnitude_sup.csv", magnitude_sup, delimiter=",")
# np.savetxt("magnitude_inf.csv", magnitude_inf, delimiter=",")

# plt.plot(magnitude_sup)
# plt.plot(magnitude_inf)
# plt.grid()
# plt.figure()

velocity_sup = get_velocity(magnitude_sup)
velocity_inf = get_velocity(magnitude_inf)
velocity_south = get_velocity(magnitude_south)
velocity_north = get_velocity(magnitude_north)
velocity_east = get_velocity(magnitude_east)
velocity_west = get_velocity(magnitude_west)

area_sup = get_area('Sup')
area_inf = get_area('Inf')
area_north = get_area('North')
area_south = get_area('South')
area_east = get_area('East')
area_west = get_area('West')

flow_sup = get_flow(velocity_sup, area_sup)
flow_inf = get_flow(velocity_inf, area_inf)
flow_north = get_flow(velocity_north, area_north)
flow_south = get_flow(velocity_south, area_south)
flow_east = get_flow(velocity_east, area_east)
flow_west = get_flow(velocity_west, area_west)

flow_in = (sum(flow_sup) + sum(flow_inf)) #+ sum(flow_west)
flow_out = (sum(flow_south) + sum(flow_east) + sum(flow_north)) #+ sum(flow_west)
mass_diff = flow_out-flow_in
print("A diferença entre os fluxos de massa é:", mass_diff)