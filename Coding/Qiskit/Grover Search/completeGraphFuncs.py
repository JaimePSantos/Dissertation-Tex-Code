import numpy as np
import matplotlib.pyplot as plt
from qiskit import *

def completeGraphWalk(N):
    qreg = QuantumRegister(N,'vertices')
    qcoin = QuantumRegister(N,'coin')
    qc = QuantumCircuit(qreg,qcoin,name='CompleteGraph')
    qc.swap(qreg[0:N],qcoin)
    return qc

def completeGraphWalk2(N):
    qreg = QuantumRegister(N,'vertices')
    qcoin = QuantumRegister(N,'coin')
    qc = QuantumCircuit(qreg,qcoin,name='CompleteGraph')
    qc.h(qcoin)
    qc.swap(qreg[0:N],qcoin)
    return qc


def hadamardCoin(N,qc,qcoin):
    qc.h(qcoin)
    return qc

def grover3Coin(N,qc,qcoin):
    qc.h(qcoin)
    qc.x(qcoin)
    qc.h(qcoin[2])
    qc.toffoli(qcoin[0],qcoin[1],qcoin[2])
    qc.h(qcoin[2])
    qc.x(qcoin)
    qc.h(qcoin)
    qc.barrier()
    return qc

def markedListComplete(markedList,N):
    oracleList = np.ones(2**N)
    for element in markedList:
        oracleList[element] = -1
    return oracleList.tolist()

def diffusionComplete(N):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(N)
    difCirc = QuantumCircuit(qreg,qcoin,name='Diffusion')
    difCirc.h(qcoin)
    
    aux = markedListComplete([0],N)
    qcAux = oracleComplete(aux,N,True)
    difCirc.append(qcAux,range(2*N))
    
    difCirc.h(qcoin)
    return difCirc

def oracleComplete(markedList,N,dif):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(N)
    qc = QuantumCircuit(qreg,qcoin,name='Oracle')
    if(dif==True):
        qc.diagonal(markedList,qcoin)
    else:
        qc.diagonal(markedList,qreg)

    return qc

def runWalkComplete(markedVertex,backend,N,times):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(N)
    creg = ClassicalRegister(N)
    qc = QuantumCircuit(qreg,qcoin,creg)
    markedVertex=markedListComplete(markedVertex,N)
    qcOracle = oracleComplete(markedVertex,N,False)
    qcDif = diffusionComplete(N)
    qcQWalk = completeGraphWalk(N)
    qc.h(qreg)

    for i in range(times):
        qc.append(qcOracle,range(2*N))
        qc.barrier()
        qc.append(qcDif,range(2*N))
        qc.barrier()
        qc.append(qcQWalk,range(2*N))
        qc.barrier()

        
    qc = transpile(qc,backend=backend,basis_gates=['cx','u3','swap'],optimization_level=3)
    qc.measure(range(N),range(N))
        
    return 
    
def runWalkComplete2(markedVertex,N,times):
    qreg = QuantumRegister(N,'vertices')
    qcoin = QuantumRegister(N,'coin')
    creg = ClassicalRegister(N)
    qc = QuantumCircuit(qreg,qcoin,creg)
    markedVertex=markedListComplete(markedVertex,N)
    qcOracle = oracleComplete(markedVertex,N,False)
    qcDif = diffusionComplete(N)
    qcQWalk = completeGraphWalk(N)
    qc.h(qreg)
    for i in range(times):
        qc.append(qcOracle,range(2*N))
        qc.barrier()
        qc.append(qcDif,range(2*N))
        qc.barrier()
        qc.append(qcQWalk,range(2*N))
        qc.barrier()

        
    qc = transpile(qc,basis_gates=['cx','u3','swap'],optimization_level=3)
    qc.measure(range(N),range(N))
        
    return qc