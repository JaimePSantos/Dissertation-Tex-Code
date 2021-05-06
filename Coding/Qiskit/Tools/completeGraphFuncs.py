import numpy as np
import matplotlib.pyplot as plt
from qiskit import *

def completeGraphWalk(N):
    qreg = QuantumRegister(N,'vertices')
    qcoin = QuantumRegister(N,'coin')
    qc = QuantumCircuit(qreg,qcoin,name='CompleteGraph')
    qc.swap(qreg[0:N],qcoin)
    qc = transpile(qc,basis_gates=['cx','u3'],optimization_level=3)
    return qc

def completeGraphWalkHCoin(N):
    qreg = QuantumRegister(N,'vertices')
    qcoin = QuantumRegister(N,'coin')
    qc = QuantumCircuit(qreg,qcoin,name='CompleteGraph')
    qc.h(qcoin)
    qc.swap(qreg[0:N],qcoin)
    return qc

def markedListComplete(markedList,N):
    oracleList = np.ones(2**N)
    for element in markedList:
        oracleList[element] = -1
    oracleList = oracleList*np.exp(1j*2*np.pi)
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

    difCirc = transpile(difCirc,basis_gates=['cx','u3'],optimization_level=3)
    return difCirc

def oracleComplete(markedList,N,dif):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(N)
    qc = QuantumCircuit(qreg,qcoin,name='Oracle')
    if(dif==True):
        qc.diagonal(markedList,qcoin)
    else:
        qc.diagonal(markedList,qreg)

    qc = transpile(qc,basis_gates=['cx','u3'],optimization_level=3)
    return qc

def runWalkComplete(markedVertex,N,backend,times):
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

        
    qc = transpile(qc,backend=backend,basis_gates=['cx','u3'],optimization_level=2)
    qc.measure(range(N),range(N))
        
    return qc
    
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
        qc.append(qcDif,range(2*N))
        qc.append(qcQWalk,range(2*N))
    qc = transpile(qc,basis_gates=['cx','u3'],optimization_level=3)
    qc.measure(range(N),range(N))
        
    return qc
