import numpy as np
import math as m
import matplotlib.pyplot as plt
from scipy import linalg
 
def AdjacencyDegreeMatrixCycle(NumberOfNodes):
    LineAdjacencyMatrix = np.zeros((NumberOfNodes,NumberOfNodes))
    LineDegreeMatrix = np.zeros((NumberOfNodes,NumberOfNodes))
    for Position in range(0,NumberOfNodes):
        LineAdjacencyMatrix[Position,(Position+1)%(NumberOfNodes)] = 1
        LineAdjacencyMatrix[(Position+1)%(NumberOfNodes),Position] = 1
        LineAdjacencyMatrix[(Position-1)%(NumberOfNodes),Position] = 1
        LineAdjacencyMatrix[Position,(Position-1)%(NumberOfNodes)] = 1
        LineDegreeMatrix[Position,Position] = 2
         
    return LineAdjacencyMatrix, LineDegreeMatrix
 
def EvolutionOperator(H,t,gamma):
    U = linalg.expm(-1j*gamma*H*t)
     
    return U
NumberOfNodes = 10
Time = 100
gamma = 1/(2*np.sqrt(2))
InitialState = np.zeros((NumberOfNodes,1))
InitialState[int(NumberOfNodes/2)] = 1
 
A, D = AdjacencyDegreeMatrixCycle(NumberOfNodes)
print(A)
U = EvolutionOperator(A,Time,gamma)
 
FinalState = U.dot(InitialState)
plt.plot(np.linspace(1,NumberOfNodes,NumberOfNodes),FinalState*np.conjugate(FinalState))
plt.show()