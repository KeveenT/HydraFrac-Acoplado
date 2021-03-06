import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import Solver
import numpy as np
from scipy import sparse
import scipy.sparse.linalg
import time
import matplotlib.pyplot as plt

class HeatTransferSolver(Solver):
    def __init__(self, workspaceDirectory, **kwargs):
        # kwargs -> outputFileName, extension, transient, verbosity
        Solver.__init__(self, workspaceDirectory, **kwargs)

    def init(self):
        self.pressureField = np.repeat(0.0, self.grid.vertices.size)
        self.prevPressureField = self.problemData.initialValues["pressure"]
        
        self.saver.save("pressure", self.prevPressureField, self.currentTime)

        self.coords,self.matrixVals = [], []
        self.difference=0.0

    def mainloop(self):
        # import leakoff
        self.assembleMatrix()
        pressureField = []
        currentTime = []
        while not self.converged:
            currentTime.append(self.currentTime)
            self.addToIndependentVector()
            self.solveLinearSystem()
            pressureField.append(self.pressureField)
            self.printIterationData()
            self.currentTime += self.timeStep
            self.saveIterationResults()
            self.checkConvergence()
            self.prevPressureField = self.pressureField
            print(currentTime[-1])
            self.iteration += 1

    def add(self, i, j, val):
        self.coords.append((i,j))
        self.matrixVals.append(val)

    def assembleMatrix(self):
        def diffusionTerm():
            # Diffusion Term
            for region in self.grid.regions:
                permeability = self.propertyData[region.handle]["Permeability"]
                viscosity = self.propertyData[region.handle]["Viscosity"]
                for element in region.elements:
                    for innerFace in element.innerFaces:
                        diffusiveFlux = (permeability/viscosity) * np.matmul( np.transpose(innerFace.globalDerivatives) , innerFace.area.getCoordinates()[:self.dimension] )
                        backwardVertexHandle = element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][0]].handle
                        forwardVertexHandle = element.vertices[element.shape.innerFaceNeighborVertices[innerFace.local][1]].handle

                        i=0
                        for vertex in element.vertices:
                            coefficient = -1.0 * diffusiveFlux[i]
                            self.add(backwardVertexHandle, vertex.handle, coefficient)
                            self.add(forwardVertexHandle, vertex.handle, -coefficient)
                            i+=1

        def accumulationTerm():
            # Transient Term
            for region in self.grid.regions:
                porosity = self.propertyData[region.handle]["Porosity"]
                fluidCompressibility = self.propertyData[region.handle]["FluidCompressibility"]
                accumulation = porosity * fluidCompressibility / self.timeStep

                for element in region.elements:
                    local = 0
                    for vertex in element.vertices:
                        self.add(vertex.handle, vertex.handle, element.subelementVolumes[local] * accumulation)
                        local += 1

        def dirichletBoundaryCondition():
            # Dirichlet Boundary Condition
            for bCondition in self.problemData.dirichletBoundaries["pressure"]:
                for vertex in bCondition.boundary.vertices:
                    self.matrixVals, self.coords = zip( *[(val, coord) for coord, val in zip(self.coords, self.matrixVals) if coord[0] != vertex.handle] )
                    self.matrixVals, self.coords = list(self.matrixVals), list(self.coords)
                    self.add(vertex.handle, vertex.handle, 1.0)

        def invertMatrix():
            # Invert Matrix
            self.matrix = sparse.csc_matrix( (self.matrixVals, zip(*self.coords)), shape=(self.grid.vertices.size, self.grid.vertices.size) )
            self.inverseMatrix = sparse.linalg.inv( self.matrix )

        diffusionTerm()
        if self.transient:
            accumulationTerm()
        dirichletBoundaryCondition()
        invertMatrix()

    def addToIndependentVector(self):
        self.independent = np.zeros(self.grid.vertices.size)

        def generationTerm():
            # Generation Term
            for region in self.grid.regions:
                fluxGeneration = self.propertyData[region.handle]["FluxGeneration"]
                for element in region.elements:
                    local = 0
                    for vertex in element.vertices:
                        self.independent[vertex.handle] += element.subelementVolumes[local] * fluxGeneration
                        local += 1

        def accumulationTerm():
            # Transient Term
            for region in self.grid.regions:
                porosity = self.propertyData[region.handle]["Porosity"]
                fluidCompressibility = self.propertyData[region.handle]["FluidCompressibility"]
                accumulation = porosity * fluidCompressibility / self.timeStep

                for element in region.elements:
                    local = 0
                    for vertex in element.vertices:
                        self.independent[vertex.handle] += element.subelementVolumes[local] * accumulation * self.prevPressureField[vertex.handle]
                        local += 1

        def neumannBoundaryCondition():
            # Neumann Boundary Condition
            for bCondition in self.problemData.neumannBoundaries["pressure"]:
                for facet in bCondition.boundary.facets:
                    for outerFace in facet.outerFaces:
                        self.independent[outerFace.vertex.handle] += bCondition.getValue(outerFace.handle) * np.linalg.norm(outerFace.area.getCoordinates())

        def leakoff():
            import PyEFVLib
            import leakoff
            global leak_sup, leak_inf
            values_sup = np.linspace(10, 1, num=30, endpoint=True)
            values_inf = np.linspace(1, 10, num=30, endpoint=True)
            # print(values_inf)
            # leak_sup = np.repeat(values_sup, 2)
            # leak_inf = np.repeat(values_inf, 2)
            leak_sup = np.full((60,), 1)
            leak_inf = np.full((60,), 1)
            grid = PyEFVLib.read("/home/keveent/PyEFVLib/meshes/msh/2D/mesh_frac.msh")
            sup = [boundary for boundary in grid.boundaries if boundary.name == "Sup"][0]
            inf = [boundary for boundary in grid.boundaries if boundary.name == "Inf"][0]
            i = 0
            for facet in sup.facets:
                for outerFace in facet.outerFaces:
                    self.independent[outerFace.vertex.handle] += leak_sup[i] * np.linalg.norm(outerFace.area.getCoordinates())
                    i = i + 1
            j = 0
            for facet in inf.facets:
                for outerFace in facet.outerFaces:
                    self.independent[outerFace.vertex.handle] += leak_inf[j] * np.linalg.norm(outerFace.area.getCoordinates())
                    j = j + 1

        def dirichletBoundaryCondition():
            # Dirichlet Boundary Condition
            for bCondition in self.problemData.dirichletBoundaries["pressure"]:
                for vertex in bCondition.boundary.vertices:
                    self.independent[vertex.handle] = bCondition.getValue(vertex.handle)

        generationTerm()
        if self.transient:
            accumulationTerm()
        neumannBoundaryCondition()
        # leakoff()
        dirichletBoundaryCondition()
        

    def solveLinearSystem(self):
        self.pressureField = np.matmul(self.inverseMatrix.toarray(), self.independent)

    def saveIterationResults(self):
        self.saver.save("pressure", self.pressureField, self.currentTime)

    def checkConvergence(self):
        self.converged = False
        # self.difference = max([abs(temp-oldTemp) for temp, oldTemp in zip(self.pressureField, self.prevPressureField)])
        self.difference = max(abs((self.pressureField - self.prevPressureField)/(max(self.pressureField) - min(self.pressureField))))
        print('Diferen??a P:', self.difference)
        if self.problemData.finalTime != None and self.currentTime > self.problemData.finalTime:
            self.converged = True
            return
        if self.iteration > 0 and self.tolerance != None and self.difference < self.tolerance:
            self.converged = True
            return
        if self.problemData.maxNumberOfIterations and self.iteration >= self.problemData.maxNumberOfIterations:
            self.converged = True
            return

def heatTransfer(workspaceDirectory, solve=True, extension="csv", saverType="default", transient=True, verbosity=True):
    solver = HeatTransferSolver(workspaceDirectory, outputFileName="Results", extension=extension, saverType=saverType, transient=transient, verbosity=verbosity)
    if solve:
        solver.solve()
    return solver

if __name__ == "__main__":
    model = "workspace/porous_flow/linear"
    if len(sys.argv)>1 and not "-" in sys.argv[1]: model=sys.argv[1]
    extension = "csv" if not [1 for arg in sys.argv if "--extension" in arg] else [arg.split('=')[1] for arg in sys.argv if "--extension" in arg][0]
    saverType = "default" if not [1 for arg in sys.argv if "--saver" in arg] else [arg.split('=')[1] for arg in sys.argv if "--saver" in arg][0]

    heatTransfer(model, extension=extension, saverType=saverType, transient=not "-p" in sys.argv, verbosity=False)
    