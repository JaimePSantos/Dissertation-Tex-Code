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
        Aer,
        transpile
        )
from qiskit.visualization import( plot_histogram,
        plot_state_city)
plt.rcParams['figure.figsize'] = 11,8
matplotlib.rcParams.update({'font.size' : 15})

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

def grover(N,steps,marked):
    qc = QuantumCircuit(N,N)
    qcOracle = oracleGrover(markedListGrover(marked,N),N)
    qcDiffusion = diffusionGrover(N)
    qc.h(range(N))
    for i in range(steps):
        qc.append(qcOracle,range(N))
        qc.barrier()
        qc.append(qcDiffusion,range(N))
        qc.barrier()
    qc.barrier()
    qc.measure(range(N),range(N))
    qc = transpile(qc,optimization_level=1)
    return qc

def saveGroverSearchFig(N,steps,markedVertex,fig, filePath, defaultFileName):
    specificFileName = ""
    i=0
    for n,m in zip(N,markedVertex):
        specificFileName+= "N%s_M%s_S"%(n,m)
        for step in steps:
            specificFileName+="%s"%step
        i+=1
        if(len(N)-i==0):
            break
        specificFileName+="_"
    savefig(fig, filePath,defaultFileName+specificFileName)
    return specificFileName

def runMultipleSearchComplete(N,steps,markedVertex):
     "Creates several instances of the coined quantum walk search circuit."
     circList = []
     circListAux = []
     for n in N:
         qreg = QuantumRegister(n)
         qsub = QuantumRegister(1)
         creg = ClassicalRegister(n)
         for step in steps:
             circ = QuantumCircuit(qreg,qsub,creg)
             circ = grover(n,step,markedVertex)
             circListAux.append(circ)
         circList.append(circListAux)
         circListAux = []
     return circList

filePath = 'GroverQiskit/'
defaultFileName = "GroverQiskitSearch_"

N = [3]
markedList = [1,2] 
steps = [0,1,2]
shots = 3000

multipleGrover = runMultipleSearchComplete(N,steps,markedList)
fig = plotMultipleQiskit(N,multipleGrover,steps,shots,True)
saveGroverSearchFig(N,steps,markedList,fig,filePath,defaultFileName)
