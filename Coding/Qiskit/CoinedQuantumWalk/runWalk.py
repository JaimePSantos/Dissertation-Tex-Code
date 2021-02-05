import sys
sys.path.append('../Tools')
from IBMTools import( 
        simul,
        savefig,
        saveMultipleHist,
        printDict,
        plotMultipleQiskit)
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from qiskit import( ClassicalRegister,
        QuantumRegister,
        QuantumCircuit,
        execute,
        Aer)
from qiskit.visualization import( plot_histogram,
                        plot_state_city)
from math import (log,ceil)
plt.rcParams['figure.figsize'] = 11,8
matplotlib.rcParams.update({'font.size' : 15})

#CNot decomposition
def cnx(qc,*qubits):
    if len(qubits) >= 3:
        last = qubits[-1]
        #A matrix: (made up of a  and Y rotation, lemma4.3)
        qc.crz(np.pi/2, qubits[-2], qubits[-1])
        #cry
        qc.cu(np.pi/2, 0, 0,0, qubits[-2],qubits[-1])
        #Control not gate
        cnx(qc,*qubits[:-2],qubits[-1])
        #B matrix (cry again, but opposite angle)
        qc.cu(-np.pi/2, 0, 0,0, qubits[-2], qubits[-1])
        #Control
        cnx(qc,*qubits[:-2],qubits[-1])
        #C matrix (final rotation)
        qc.crz(-np.pi/2,qubits[-2],qubits[-1])
    elif len(qubits)==3:
        qc.ccx(*qubits)
    elif len(qubits)==2:
        qc.cx(*qubits)
    return qc

def incr(qwc,q,subnode,n):
    for j in range(-1,n-1):
        if(j==-1):
            cnx(qwc,subnode[0],*q[-1::-1])
            qwc.barrier()
        else:
            cnx(qwc,subnode[0],*q[-1:j:-1])
            qwc.barrier()
    return qwc

def decr(qwc,q,subnode,n):
    qwc.x(subnode[0])
    c=0
    qwc.x(q[-1:0:-1])
    for j in range(-1,n-1):
        if(j==-1):
            c+=1
            cnx(qwc,subnode[0],*q[-1::-1])
            qwc.x(q[c])
            qwc.barrier()
        else:
            c+=1
            cnx(qwc,subnode[0],*q[-1:j:-1])
            if(c==n):
                break
            qwc.x(q[c])
            qwc.barrier()
    qwc.x(subnode[0])
    return qwc

def runWalk(N,steps,stateVec):
    "Creates a single instance of the coined quantum walk cicuit."
    qreg = QuantumRegister(N)
    qsub = QuantumRegister(1)
    creg = ClassicalRegister(N)
    qwc = QuantumCircuit(qreg,qsub,creg)
    qwc.x(qreg[0])
    for i in range(0,steps):
        qwc.h(qsub[0])
        incr(qwc,qreg,qsub,N)
        decr(qwc,qreg,qsub,N)
    if not stateVec:
        qwc.barrier()
        qwc.measure(qreg,creg)
        qwc.barrier()
    return qwc

def saveCoinedWalkFig(N,steps,fig, filePath, defaultFileName):
    specificFileName = ""
    i=0
    for n in N:
        specificFileName+= "N%s_S"%n
        for step in steps:
            specificFileName+="%s"%step
        i+=1
        if(len(N)-i==0):
            break
        specificFileName+="_"
    savefig(fig, filePath,defaultFileName+specificFileName)
    return specificFileName

def runMultipleWalks(N,steps,stateVec):
    "Creates several instances of the coined quantum walk circuit."
    circList = []
    circListAux = []
    for n in N:
        qreg = QuantumRegister(n)
        qsub = QuantumRegister(1)
        creg = ClassicalRegister(n)
        for step in steps:
            circ = QuantumCircuit(qreg,qsub,creg)
            circ = runWalk(n,step,stateVec)
            circListAux.append(circ)
        circList.append(circListAux)
        circListAux = []
    return circList



filePath = 'CoinedQuantumWalk/'
defaultFileName = "CoinedQW_"

singleN = 3 
singleSteps = 1 

N=[4]
steps=[0,1,2,3,4,5]
shots = 3000
multipleWalks = runMultipleWalks(N,steps,False)
fig = plotMultipleQiskit(N,multipleWalks,steps,shots,True)
saveCoinedWalkFig(N,steps,fig,filePath,defaultFileName)
