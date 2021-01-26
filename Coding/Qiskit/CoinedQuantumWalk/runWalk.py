import sys
sys.path.append('../Tools')
from IBMTools import( 
        simul,
        savefig,
        saveMultipleHist,
        printDict)
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

def runWalk(N,times,stateVec):
    qreg = QuantumRegister(N)
    qsub = QuantumRegister(1)
    creg = ClassicalRegister(N)
    qwc = QuantumCircuit(qreg,qsub,creg)
    for i in range(0,times):
        qwc.h(qsub[0])
        incr(qwc,qreg,qsub,N)
        decr(qwc,qreg,qsub,N)
    if not stateVec:
        qwc.barrier()
        qwc.measure(qreg,creg)
        qwc.barrier()
    return qwc

def baseResultDict(n):
    baseDict = {}
    for decNumber in range(2**n):
        decToBin = bin(decNumber)[2:].zfill(ceil(log(2**n,2)))
        baseDict[str(decToBin)] = 0
    return baseDict

def multBaseResultDict(N,steps):
    baseResultDictList = []
    for n in N:
        for step in steps:
            baseDict = baseResultDict(n)
            baseResultDictList.append(baseDict)
    return baseResultDictList

def multNormalizedResultDict(baseDictList,qiskitDictList):
    normalizedResultDictList = []
    for baseDict,qiskitDict in zip(baseDictList,qiskitDictList):
        normalizedResultDict = {**baseDict,**qiskitDict}
        normalizedResultDictList.append(normalizedResultDict)
    return normalizedResultDictList
#TODO: Fazer merge dos dois dicionarios.
def multResultsSim(multipleCircs,shots):
    resultList = []
    result = {}
    correctedResult = {}
    for circList in multipleCircs:
        for circ in circList:
            result = simul(circ,False,shots)
            correctedResult = { k[::-1] : v/shots for k, v in result.items()}
            resultList.append(correctedResult)
            result = {}
    return resultList

#TODO: Decidir como ajustar os eixos da figura. 
#TODO: Decidir como por os labels dos eixos da figura.
def multSubPlot(resultList):
    nrows = len(resultList) 
    ncols = 1
    index = 1
    fig = plt.figure()
    axs = []
    for resultAux in resultList:
        axs.append(fig.add_subplot(nrows,ncols,index))
        axs[-1].bar(resultAux.keys(),resultAux.values(),width=0.4)
        index+=1
    for ax in axs:
        axs[-1].get_shared_y_axes().join(axs[-1],ax)
    for ax in axs[:-1]:
        ax.set_xticklabels([])
    fig.tight_layout(pad=1.0)
    return fig

def plotMultipleQiskit(N,multipleCircs,steps,shots):
    qiskitResultList = multResultsSim(multipleCircs,shots)
    baseDictList = multBaseResultDict(N,steps)
    normalizedResultDictList = multNormalizedResultDict(baseDictList,qiskitResultList)
    fig = multSubPlot(normalizedResultDictList)
    plt.show()

def runMultipleWalks(N,steps,stateVec):
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
defaultFileName = "CoinedQW_N"

singleN = 3 
singleSteps = 1 

N=[3]
steps=[0,1,2]
shots = 3000
multipleWalks = runMultipleWalks(N,steps,False)
plotMultipleQiskit(N,multipleWalks,steps,shots)
