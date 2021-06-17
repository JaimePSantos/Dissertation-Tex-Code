import numpy as np
import matplotlib as mpl 
mpl.use('TkAgg')
import matplotlib.pyplot as plt
from qiskit import( ClassicalRegister,
        QuantumRegister,
        QuantumCircuit,
        execute,
        Aer,
        IBMQ,
        transpile)
from qiskit.tools.monitor import job_monitor
from qiskit.providers.ibmq import least_busy
from qiskit.providers.aer.noise import NoiseModel
from qiskit.visualization import( plot_histogram,
                        plot_state_city,
                        plot_gate_map, 
                        plot_circuit_layout,
                        circuit_drawer)
from qiskit.circuit.library import QFT
from math import (log,ceil)
from scipy.fft import fft, ifft
from scipy.linalg import dft, inv, expm, norm
from numpy.linalg import matrix_power
#import networkx as nx
mpl.rcParams['figure.figsize'] = 11,8
mpl.rcParams.update({'font.size' : 15})

def simul(qc,stateVec,shots):
    if stateVec:
        backend = Aer.get_backend('statevector_simulator')
        result = execute(qc,backend,shots=shots).result().get_statevector(qc,decimals=3)
    else:
        backend = Aer.get_backend('qasm_simulator')
        result = execute(qc,backend,shots=shots).result().get_counts()
    return result

def decResultDict(n):
    "Retuns a dictionary composed of a range of N keys converted to binary."
    baseDict = {}
    for decNumber in range(2**n):
        dec = decNumber 
        baseDict[dec] = 0
    return baseDict

def normalizedResults(resultsDict,n,shots):
    decDict = decResultDict(n)
    correctedResults = {int(k,2) : v/shots for k,v in resultsDict.items()}
    newDict1 = correctedResults
    newDict2 = decDict
    normalizedResults = {**newDict2,**newDict1}
    return normalizedResults

def c_increment(n):
    c_inc = QuantumCircuit(n)
    controls = [x for x in range(n-1)]
    for p in range(n-1):
        c_inc.mcx(controls,controls[-1] + 1)
        controls.pop()
    c_inc.x(0)
    return c_inc
    
def c_decrement(n):
    c_dec = QuantumCircuit(n)
    controls = [x for x in range(n-1)]
    c_dec.x(controls)
    for p in range(n-2):
        c_dec.mcx(controls,controls[-1] + 1)
        c_dec.x(controls[-1])
        controls.pop()   
    c_dec.cx(0,1)
    return c_dec

def initialCond(qc,string,N):
    for x in range(N):
        if string[x] == '1':
            qc.x(x)
    return qc

def stagWalk(N,theta,steps,initString):
    qreg = QuantumRegister(N)
    creg = ClassicalRegister(N)
    qc = QuantumCircuit(qreg,creg)
    qc = initialCond(qc,initString,N)
    qcInc = c_increment(N)
    qcDec = c_decrement(N)
    for step in range(steps):
        qc.rx(2*theta,qreg[0])
        qc.barrier()
        qc.append(qcInc,qreg)
        qc.barrier()
        qc.rx(2*theta,qreg[0])
        qc.barrier()
        qc.append(qcDec,qreg)
        qc.barrier()
    qc.measure(qreg,creg)
    return qc

def multContCirc(N,theta,stepList,initString):
    circList = []
    for steps in zip(stepList):
        circ =  stagWalk(N,theta,steps,initString)
        #circ = transpile(circ,optimization_level=3,backend=backend, layout_method=method)
        circList.append(circ)
    return circList

N=3
theta = np.pi/3
steps = 1
initString = '001'
stagQC = stagWalk(N,theta,steps,initString)
stagQC = stagQC.decompose()
stagQC.draw(output='mpl')
plt.show() 