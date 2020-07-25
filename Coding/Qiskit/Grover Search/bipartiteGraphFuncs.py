import numpy as np
import matplotlib.pyplot as plt
from qiskit import *

def bipartiteWalk(N,qc,qreg,qcoin):
    qc.x(qreg[N-1])
    qc.swap(qreg[:-1],qcoin)
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

def newWalk(N,qc,qreg,qcoin):
    qc.swap(qreg[0],qcoin[0])
    qc.swap(qreg[1],qcoin[1])    
    qc.swap(qreg[2],qcoin[2])
    qc.x(qreg[3])
    return qc

def oracle(marked,N):
    qreg = QuantumRegister(N)
    qc = QuantumCircuit(qreg)
    D= np.ones(2**N)
    D[marked] = -1
    #D[63] = -1
    #D[62] = -1
    D[1] = -1
    D[15]=-1
    #D[61] = -1
    #D[60] = -1
    D[14] = -1

    print(D)
    qc.diagonal(D.tolist(),qreg)
    qc = transpile(qc,basis_gates =['cx','u3'],optimization_level=3)
    return qc

def runWalk(qc,qreg,qcoin,markedVertex,Ncoin,N,times):
    qcaux = oracle(markedVertex,N)
    qc.h(qreg)
    
    for i in range(times):
        qc.append(qcaux,range(N))
        qc.barrier()
        qc= grover3Coin(Ncoin,qc,qcoin)
        qc= newWalk(N,qc,qreg,qcoin)
        qc.barrier()
        
    return qc