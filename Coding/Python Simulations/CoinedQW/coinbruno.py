import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import math as m
 
def CoinOperator(CoinString):
    if CoinString == "X":
        Coin = np.array([[0,1],[1,0]])
    elif CoinString == "H":
        Coin = np.array([[1,1],[1,-1]])/np.sqrt(2)
         
    return Coin
     
def CycleShiftOperator(NumberOfNodes):
    ShiftPlus = np.zeros((NumberOfNodes,NumberOfNodes))
    ShiftMinus = np.zeros((NumberOfNodes,NumberOfNodes))
     
    for Node in range(0,NumberOfNodes):
        ShiftPlus[(Node+1)%(NumberOfNodes),Node] = 1
        ShiftMinus[(Node-1)%(NumberOfNodes),Node] = 1
         
    ShiftOperator = np.kron(np.array([[1,0],[0,0]]),ShiftPlus) + np.kron(np.array([[0,0],[0,1]]),ShiftMinus)
    print("Shift matrix:\n",ShiftOperator)
     
    return ShiftOperator
 
def CoinedUnitaryOperator(Coin,Shift):
    UnitaryOperator = Shift.dot(np.kron(Coin,np.eye(NumberOfNodes)))
     
    return UnitaryOperator
 
def GraphInitialState(NumberOfNodes,Positions,Amplitudes,CoinState):
    InitialState = np.zeros((NumberOfNodes,1))
    for Node in range(0,Positions.size):
        InitialState[Positions[Node],0] = Amplitudes[Node]
    print("initstate")
    print(InitialState)
         
    return np.kron(CoinState,InitialState)
 
def Evolution(U,InitialState,Steps):
    for t in range(0,Steps):
        InitialState = U.dot(InitialState)
         
    return InitialState
 
def CycleProbability(FinalState,NumberOfNodes):
    ProbabilityVector = np.zeros((NumberOfNodes,1))
    for x in range(0,NumberOfNodes):
        ProbabilityVector[x] = FinalState[x]*np.conjugate(FinalState[x]) + FinalState[NumberOfNodes + x]*np.conjugate(FinalState[NumberOfNodes + x])
 
    return ProbabilityVector
     
NumberOfNodes = 6
Steps = 3
 
Coin = CoinOperator("H")
# print(Coin)
# print("##################################################")
Shift = CycleShiftOperator(NumberOfNodes)
# print(Shift)
# print("##################################################")
U = CoinedUnitaryOperator(Coin,Shift)
# print(U)
# print("UUUUUUUU##################################################")
 
Amplitudes = np.array([1])
blabla = np.array([[1/np.sqrt(2)],[(1*1j)/np.sqrt(2)]])
# print(blabla)
# print("coinstate##################################################")
InitialState = GraphInitialState(NumberOfNodes,np.array([int((NumberOfNodes+1)/2)]),Amplitudes,np.array([[1/np.sqrt(2)],[(1*1j)/np.sqrt(2)]]))
# print(InitialState)
# print("psi0##################################################")
FinalState = Evolution(U,InitialState,Steps)
# print(FinalState)
# print("##################################################")

ProbabilityVector = CycleProbability(FinalState,NumberOfNodes)
 
plt.plot(np.linspace(1,NumberOfNodes,NumberOfNodes),ProbabilityVector)
#plt.show()