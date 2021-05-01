import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__), os.path.pardir))
from PyEFVLib import Solver
import numpy as np
from scipy import sparse
import scipy.sparse.linalg
import time
from PyEFVLib import MeshioSaver
import mesh_reader
import frac_flow_NS
import leakoff
import csv
import datetime

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
        global pressureField, Sfrac, S_Old, S, Pf, Pf_Old, currentTime, leak_sup, leak_inf, leak_sup_old, leak_inf_old, Pn, Ps
        begin_time = datetime.datetime.now()
        pressureField = []
        pressureField.append(np.array(self.prevPressureField, dtype=float))

        #Navier-Stokes
        nx = mesh_reader.get_nx() #Nบmero de Volumes
        nu = (nx + 1) #Nรบmero de volumes de velocidade
        S = []
        Pf = []
        Sfrac = []
        S.append(frac_flow_NS.get_initial_S())
        Sfrac.append(frac_flow_NS.get_initial_S())
        S_Old = S[-1]
        Pn, Ps = leakoff.get_Pns(pressureField[-1])
        S.append(frac_flow_NS.solve(S_Old, Pn, Ps))
        Pf.append(S[-1][nu:])

        leak_sup_old, leak_inf_old = leakoff.get_leakoff_frac(Pf, Pn, Ps)
        leakoff_iterator = []
        self.assembleMatrix()
        checkS = 'Not Converged'
        while (not self.converged):
            print('Tempo atual:', self.currentTime)
            check = 'Not Converged'
            iterator = 0
            while check != 'Converged':
                self.addToIndependentVector()
                self.solveLinearSystem()
                pressureField[-1] = self.pressureField
                Pn, Ps = leakoff.get_Pns(pressureField[-1])
                leak_sup, leak_inf = leakoff.get_leakoff_frac(Pf, Pn, Ps)
                check = leakoff.leakoff_check(iterator, leak_sup, leak_inf, leak_sup_old, leak_inf_old)
                leak_sup_old = leak_sup.copy()
                leak_inf_old = leak_inf.copy()
                S[-1] = frac_flow_NS.solve(S_Old, Pn, Ps) #NS
                Pf.append(S[-1][nu:]) #NS
                iterator = iterator + 1
            leakoff_iterator.append(iterator)
            self.currentTime += self.timeStep
            self.saveIterationResults()
            self.checkConvergence()
            checkS = frac_flow_NS.checkS(S_Old, S) #NS
            Sfrac.append(S[-1]) #NS
            S_Old = S[-1] #NS
            self.prevPressureField = self.pressureField
            self.iteration += 1
        print(datetime.datetime.now() - begin_time)
        # runtime = datetime.datetime.now() - begin_time
        # np.savetxt("runtime.txt", runtime, delimiter=",")
        np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/no_leak_iterations_ns.csv", leakoff_iterator, delimiter=",")
        np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/leakoff_superior_ns.csv", leak_sup, delimiter=",")
        np.savetxt("/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/leakoff_inferior_ns.csv", leak_inf, delimiter=",")
        with open('/home/keveent/PyEFVLib/Hydraulic_Fracture/Resultados/navier_stokes.csv', 'w') as csv_file:
            csv_writer = csv.writer(csv_file)
            csv_writer.writerows(Sfrac)
        return

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

        def leakoff(): #Adiciona o leakoff ao termo fonte como uma condição de Neumann
            sup = [boundary for boundary in self.grid.boundaries if boundary.name == "Sup"][0]
            inf = [boundary for boundary in self.grid.boundaries if boundary.name == "Inf"][0]
            i = 0
            for facet in sup.facets:
                for outerFace in facet.outerFaces:
                    self.independent[outerFace.vertex.handle] += leak_sup_old[i] * np.linalg.norm(outerFace.area.getCoordinates())
                    i = i + 1
            j = 0
            for facet in inf.facets:
                for outerFace in facet.outerFaces:
                    self.independent[outerFace.vertex.handle] += leak_inf_old[j] * np.linalg.norm(outerFace.area.getCoordinates())
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
        leakoff()
        dirichletBoundaryCondition()
        
    def solveLinearSystem(self):
        self.pressureField = np.matmul(self.inverseMatrix.toarray(), self.independent)

    def saveIterationResults(self):
        self.saver.save("pressure", self.pressureField, self.currentTime)

    def checkConvergence(self):
        self.converged = False
        # self.difference = max([abs((temp-oldTemp)) for temp, oldTemp in zip(self.pressureField, self.prevPressureField)])
        self.difference = max(abs((self.pressureField - self.prevPressureField)/(max(self.pressureField) - min(self.pressureField))))
        print('Diferença P:', self.difference)
        if self.problemData.finalTime != None and self.currentTime > self.problemData.finalTime:
            self.converged = True
            return
        if self.iteration > 0 and self.tolerance != None and self.difference < self.tolerance:
            self.converged = True
            print("Pressão no Meio Convergida")
            return
        if self.problemData.maxNumberOfIterations and self.iteration >= self.problemData.maxNumberOfIterations:
            self.converged = True
            return

def heatTransfer(workspaceDirectory, solve=True, extension="csv", saverType="default", transient=True, verbosity=True):
    solver = HeatTransferSolver(workspaceDirectory, outputFileName="Results_ns", extension=extension, saverType=saverType, transient=transient, verbosity=verbosity)
    if solve:
        solver.solve()
    return solver

if __name__ == "__main__":
    model = "workspace/porous_flow/linear"
    if len(sys.argv)>1 and not "-" in sys.argv[1]: model=sys.argv[1]
    extension = "csv" if not [1 for arg in sys.argv if "--extension" in arg] else [arg.split('=')[1] for arg in sys.argv if "--extension" in arg][0]
    saverType = "default" if not [1 for arg in sys.argv if "--saver" in arg] else [arg.split('=')[1] for arg in sys.argv if "--saver" in arg][0]

    heatTransfer(model, extension=extension, saverType=saverType, transient=not "-p" in sys.argv, verbosity=False)
    #heatTransfer(model, extension="xdmf", saverType=saverType, transient=not "-p" in sys.argv, verbosity=False)

