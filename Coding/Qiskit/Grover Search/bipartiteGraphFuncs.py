import numpy as np
import matplotlib.pyplot as plt
from qiskit import *

def bipartiteWalk(N,n,qreg,qcoin):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(n)
    qc = QuantumCircuit(qreg,qcoin,name='BipartiteGraph')
    qc.x(qreg[N-1])
    qc.swap(qreg[0:N-1],qcoin[0:n])
    return qc
    
def markedListBipartite(markedList,N,n,dif):
    if(dif==False):    
        oracleList = np.ones(2**N)
        for element in markedList:
            oracleList[element] = -1
    if(dif==True):
        oracleList = np.ones(2**n)
        for element in markedList:
            oracleList[element] = -1

    return oracleList.tolist()

def diffusionBipartite(N,n):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(n)
    difCirc = QuantumCircuit(qreg,qcoin,name='Diffusion')
    difCirc.h(qcoin)
    
    aux = markedListBipartite([0],N,n,True)
    qcAux = oracleBipartite(aux,N,n,True)
    difCirc.append(qcAux,range(n+N))
    difCirc.h(qcoin)

    difCirc = transpile(difCirc,basis_gates=['cx','H','swap','u3','x'],optimization_level=3)
    return difCirc

def oracleBipartite(markedList,N,n,dif):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(n)
    qc = QuantumCircuit(qreg,qcoin,name='Oracle')
    if(dif==True):
        qc.diagonal(markedList,qcoin)
    else:
        qc.diagonal(markedList,qreg)
    qc = transpile(qc,basis_gates=['cx','H','swap','u3','x'],optimization_level=3)
    return qc

def runWalkBipartite2(markedVertex,N,n,times):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(n)
    creg = ClassicalRegister(N)
    qc = QuantumCircuit(qreg,qcoin,creg)
    markedVertex=markedListBipartite(markedVertex,N,n,False)
    qcOracle = oracleBipartite(markedVertex,N,n,False)
    qcDif = diffusionBipartite(N,n)
    qcQWalk = bipartiteWalk(N,n,qreg,qcoin)
    qc.h(qreg)
    for i in range(times):
        qc.append(qcOracle,range(n+N))
        qc.barrier()
        qc.append(qcDif,range(n+N))
        qc.barrier()
        qc.append(qcQWalk,range(n+N))
        qc.barrier()

        
    qc = transpile(qc,basis_gates=['cx','H','swap','u3','x'],optimization_level=3)
    qc.measure(range(N),range(N))
        
    return qc