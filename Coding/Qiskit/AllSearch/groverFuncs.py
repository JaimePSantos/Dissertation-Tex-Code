import numpy as np
import matplotlib.pyplot as plt
from qiskit import *
from qiskit.visualization import plot_histogram
from qiskit import IBMQ
from qiskit.tools.monitor import job_monitor

def markedListGrover(markedList,N):
    oracleList = np.ones(2**N)
    for element in markedList:
        oracleList[element] = -1
    return oracleList.tolist()

def getOracle(markedList,N):
    oracleList = np.eye(2**N)
    for element in markedList:
        oracleList[element][element] = -1
    return oracleList

def oracleGrover(markedList,N):
    qreg = QuantumRegister(N)
    qc = QuantumCircuit(qreg,name='Oracle')
    qc.diagonal(markedList,qreg)

    qc=transpile(qc,optimization_level=3)
    return qc


def diffusionGrover(N):
    qreg = QuantumRegister(N)
    difCirc = QuantumCircuit(qreg,name='Diffusion')
    difCirc.h(qreg)
    aux = markedListGrover([0],N)
    qcAux = oracleGrover(aux,N)
    difCirc.append(qcAux,range(N))
    difCirc.h(qreg)
    difCirc=transpile(difCirc,optimization_level=3)

    
    return difCirc

def grover(marked,N,backend,steps):
    qc = QuantumCircuit(N,N)
    qcOracle = oracleGrover(markedListGrover(marked,N),N)
    qcDiffusion = diffusionGrover(N)
    qc.h(range(N))
    for i in range(steps):
        qc.append(qcOracle,range(N))
        qc.barrier()
        qc.append(qcDiffusion,range(N))
        qc.barrier()
    qc = transpile(qc,basis_gates=['cx','u3','u2','u1','id'],backend=backend,optimization_level=3)
    qc.barrier()
    qc.measure(range(N),range(N))
    return qc

def grover2(marked,N,steps):
    qc = QuantumCircuit(N,N)
    qcOracle = oracleGrover(markedListGrover(marked,N),N)
    qcDiffusion = diffusionGrover(N)
    qc.h(range(N))
    for i in range(steps):
        qc.append(qcOracle,range(N))
        qc.barrier()
        qc.append(qcDiffusion,range(N))
        qc.barrier()
    #qc = transpile(qc,basis_gates=['cx','u3'],optimization_level=3)
    qc.barrier()
    qc.measure(range(N),range(N))
    return qc

def simul(qc):
    backend = Aer.get_backend('qasm_simulator')
    result = execute(qc,backend,shots=3000).result().get_counts()
    return result

def simulUnitary(qc):
    backend = Aer.get_backend('unitary_simulator')
    result = execute(qc,backend).result()
    unitary = result.get_unitary(qc)
    return unitary