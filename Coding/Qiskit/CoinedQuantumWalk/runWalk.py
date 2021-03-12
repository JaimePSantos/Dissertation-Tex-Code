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
    if len(qubits) > 3:
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
            #qwc.barrier()
        else:
            cnx(qwc,subnode[0],*q[-1:j:-1])
           # qwc.barrier()
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
            #qwc.barrier()
        else:
            c+=1
            cnx(qwc,subnode[0],*q[-1:j:-1])
            if(c==n):
                break
            qwc.x(q[c])
            #qwc.barrier()
    qwc.x(subnode[0])
    return qwc

def incrCirc(qc,q,subnode,n,toGate):
    for j in range(-1,n-1):
        if(j==-1):
            cnx(qc,subnode[0],*q[-1::-1])
        else:
            cnx(qc,subnode[0],*q[-1:j:-1])
    if toGate:
        qc = qc.to_gate()
        qc.name = '      INC      '
    return qc

def decrCirc(qc,q,subnode,n,toGate):
    qc.x(subnode[0])
    c=0
    qc.x(q[-1:0:-1])
    for j in range(-1,n-1):
        if(j==-1):
            c+=1
            cnx(qc,subnode[0],*q[-1::-1])
            qc.x(q[c])
        else:
            c+=1
            cnx(qc,subnode[0],*q[-1:j:-1])
            if(c==n):
                break
            qc.x(q[c])
    qc.x(subnode[0])
    if toGate:
        qc = qc.to_gate()
        qc.name = '      DEC      '
    return qc

def runWalk(N,steps,stateVec):
    "Creates a single instance of the coined quantum walk cicuit."
    qreg = QuantumRegister(N)
    qsub = QuantumRegister(1)
    creg = ClassicalRegister(N)
    qwc = QuantumCircuit(qreg,qsub,creg)
    qwc.x(qreg[0])
    for i in range(0,steps):
        qwc.h(qsub[0])
        qwc.barrier()
        incr(qwc,qreg,qsub,N)
        qwc.barrier()
        decr(qwc,qreg,qsub,N)
        qwc.barrier()
    if not stateVec:
        qwc.measure(qreg,creg)
    return qwc

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


def circRunWalk(N,steps):
    "Creates a single instance of the coined quantum walk cicuit."
    qreg = QuantumRegister(N,name='node')
    qsub = QuantumRegister(1, name='coin')
    creg = ClassicalRegister(N)
    qwc = QuantumCircuit(qreg,qsub,creg)
    incrCirc1 = QuantumCircuit(qreg,qsub)
    decrCirc1 = QuantumCircuit(qreg,qsub)
    incrCirc1 = incrCirc(incrCirc1,qreg,qsub,N,True)
    decrCirc1 = decrCirc(decrCirc1,qreg,qsub,N,True)
    qwc.x(qreg[0])
    qwc.h(qsub[0])
    qwc.barrier()
    for i in range(0,steps):
        qwc.append(incrCirc1,[N]+list(range(0,N)))
        qwc.append(decrCirc1,[N]+list(range(0,N)))
        qwc.barrier()
        if i!=steps-1:
            qwc.h(qsub[0])
    qwc.measure(qreg,creg)
    return qwc

def circRunMultipleWalks(N,steps):
    "Creates several instances of the coined quantum walk circuit."
    circList = []
    circListAux = []
    for n in N:
        qreg = QuantumRegister(n)
        qsub = QuantumRegister(1)
        creg = ClassicalRegister(n)
        for step in steps:
            circ = QuantumCircuit(qreg,qsub,creg)
            circ = circRunWalk(n,step)
            circListAux.append(circ)
        circList.append(circListAux)
        circListAux = []
    return circList

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

def printIncr(N,steps,style):
    "Creates a single instance of the coined quantum walk cicuit."
    for n in N:
        qreg = QuantumRegister(n,name='node')
        qsub = QuantumRegister(1, name='coin')
        incrCirc1 = QuantumCircuit(qreg,qsub)
        incrCirc1 = incrCirc(incrCirc1,qreg,qsub,n,False)
        fig = incrCirc1.draw(output='mpl',style=style) 
    return fig 

def printDecr(N,steps,style):
    "Creates a single instance of the coined quantum walk cicuit."
    for n in N:
        qreg = QuantumRegister(n,name='node')
        qsub = QuantumRegister(1, name='coin')
        decrCirc1 = QuantumCircuit(qreg,qsub)
        decrCirc1 = decrCirc(decrCirc1,qreg,qsub,n,False)
        fig = decrCirc1.draw(output='mpl',style=style) 
    return fig 


def drawCirc(circMultWalk,style):
    for circList in circMultWalk:
        for circ in circList:
            fig = circ.draw(output='mpl',style=style)
    return fig


filePath = 'CoinedQuantumWalk/'
circFilePath = 'CoinedQuantumWalk/Circuits/'
circIncrFilePath = 'CoinedQuantumWalk/Circuits/'
circDecrFilePath = 'CoinedQuantumWalk/Circuits/'
defaultFileName = "CoinedQW_"
circDefaultFileName = "circCoinedQW_"
circIncrDefaultFileName = "circIncr_"
circDecrDefaultFileName = "circDecr_"
style = {'figwidth':18,'fontsize':17,'subfontsize':14}
styleIncr = {'figwidth':15,'fontsize':17,'subfontsize':14, 'compress':True}
styleDecr = {'figwidth':15,'fontsize':17,'subfontsize':14 }

#singleN = 3 
#singleSteps = 1 

#Coined quantum walk probability distribution.
N=[3]
steps=[0,1,2,3]
shots = 3000
multipleWalks = runMultipleWalks(N,steps,False)
fig = plotMultipleQiskit(N,multipleWalks,steps,shots,True)
saveCoinedWalkFig(N,steps,fig,filePath,defaultFileName)

#Coined quantum walk probability distribution.
N5=[5]
steps5=[0,5,10,15,20]
shots = 3000
multipleWalks5 = runMultipleWalks(N5,steps5,False)
fig5 = plotMultipleQiskit(N5,multipleWalks5,steps5,shots,True)
saveCoinedWalkFig(N5,steps5,fig5,filePath,defaultFileName)

#Coined quantum walk circuit.
circN = [3]
circSteps = [3]
circMultWalk = circRunMultipleWalks(circN,circSteps)
circFig = drawCirc(circMultWalk,style)
saveCoinedWalkFig(circN,circSteps,circFig,circFilePath,circDefaultFileName)

#Increment circuit.
circIncrN = [3]
circIncrSteps = [3]
circIncrFig= printIncr(circIncrN,circIncrSteps,styleIncr)
saveCoinedWalkFig(circIncrN,circIncrSteps,circIncrFig,circIncrFilePath,circIncrDefaultFileName)

#Decrement circuit.
circDecrN = [3]
circDecrSteps = [3]
circDecrFig= printDecr(circDecrN,circDecrSteps,styleDecr)
saveCoinedWalkFig(circDecrN,circDecrSteps,circDecrFig,circDecrFilePath,circDecrDefaultFileName)
