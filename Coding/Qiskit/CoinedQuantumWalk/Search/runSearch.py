import sys
sys.path.append('../../Tools')
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

def bipartiteWalk(N,n,qreg,qcoin):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(n)
    qc = QuantumCircuit(qreg,qcoin,name='BipartiteGraph')
    qc.x(qreg[N-1])
    qc.swap(qreg[0:N-1],qcoin[0:n])
    return qc

def completeGraphWalk(N):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(N)
    qc = QuantumCircuit(qreg,qcoin,name='CompleteGraph')
    qc.swap(qreg[0:N],qcoin)
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

def drawDiffusionComplete(N):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(N)
    difCirc = QuantumCircuit(qreg,qcoin,name='     Diff     ')
    difCirc.h(qcoin)
    aux = markedListComplete([0],N)
    qcAux = oracleComplete(aux,N,True)
    difCirc.append(qcAux,range(2*N))
    difCirc.h(qcoin)
    difCirc = transpile(difCirc)#,basis_gates=['cx','u3'],optimization_level=3)
    return difCirc

def oracleComplete(markedList,N,dif):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(N)
    qc = QuantumCircuit(qreg,qcoin,name='    Oracle     ')
    if(dif==True):
        qc.diagonal(markedList,qcoin)
    else:
        qc.diagonal(markedList,qreg)
    qc = transpile(qc,basis_gates=['cx','u3'],optimization_level=3)
    return qc

def drawOracleComplete(markedList,N,dif):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(N)
    qc = QuantumCircuit(qreg,qcoin,name='    Oracle     ')
    if(dif==True):
        qc.diagonal(markedList,qcoin)
    else:
        qc.diagonal(markedList,qreg)
    qc = transpile(qc)#,basis_gates=['cx','u3'],optimization_level=3)
    return qc

def runSearchComplete(N,steps,markedVertex):
    qreg = QuantumRegister(N,'vertices')
    qcoin = QuantumRegister(N,'coin')
    creg = ClassicalRegister(N)
    qc = QuantumCircuit(qreg,qcoin,creg)
    markedVertex=markedListComplete(markedVertex,N)
    qcOracle = oracleComplete(markedVertex,N,False)
    qcDif = diffusionComplete(N)
    qcQWalk = completeGraphWalk(N)
    qc.h(qreg)
    for i in range(steps):
        qc.append(qcOracle,range(2*N))
        qc.append(qcDif,range(2*N))
        qc.append(qcQWalk,range(2*N))
    qc = transpile(qc,basis_gates=['cx','u3'],optimization_level=1)
    qc.measure(range(N),range(N))
    return qc

def saveCoinedSearchFig(N,steps,markedVertex,fig, filePath, defaultFileName):
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
            circ = runSearchComplete(n,step,markedVertex)
            circListAux.append(circ)
        circList.append(circListAux)
        circListAux = []
    return circList

def drawSearchComplete(N,steps,markedVertex,style):
    qreg = QuantumRegister(N,'qv')
    qcoin = QuantumRegister(N,'qc')
    creg = ClassicalRegister(N)
    qc = QuantumCircuit(qreg,qcoin,creg)
    markedVertex=markedListComplete(markedVertex,N)
    qcOracle = drawOracleComplete(markedVertex,N,False)
    qcDif = drawDiffusionComplete(N)
    qcQWalk = completeGraphWalk(N)
    qc.h(qreg)
    for i in range(steps):
        qc.append(qcOracle,range(2*N))
        qc.append(qcDif,range(2*N))
        qc.append(qcQWalk,range(2*N))
    qc.measure(range(N),range(N))
    qc = transpile(qc)
    fig = qc.draw(output='mpl',style=style)
    return fig

filePath = 'CoinedQuantumWalk/Search/'
defaultFileName = "CoinedQiskitSearch_"
circFilePath = 'CoinedQuantumWalk/Search/Circuits'
defaultCircFileName = "GroverQiskitCirc_"
defaultCircOracleFileName = "GroverQiskitCircOracle_"
defaultCircDiffFileName = "GroverQiskitCircDiff_"


style = {'figwidth':20,'fontsize':17,'subfontsize':14}#,'compress':True}
 
drawSearchComplete(3,3,[0],style)
plt.show()
##TODO: markedVertex labels are not correct due to post processing.
#N=[4]
#steps=[0,1,2,3,4]
#markedVertex = [1] 
#shots = 3000
#multipleWalks = runMultipleSearchComplete(N,steps,markedVertex)
#fig = plotMultipleQiskit(N,multipleWalks,steps,shots,True)
#saveCoinedSearchFig(N,steps,markedVertex,fig,filePath,defaultFileName)
