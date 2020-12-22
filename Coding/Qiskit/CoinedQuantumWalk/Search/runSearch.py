import sys
sys.path.append('../../Tools')
from IBMTools import( 
        simul,
        savefig)
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
        qc.measure(range(N),range(N))
    qc = transpile(qc,basis_gates=['cx','u3'],optimization_level=1)
    qc.measure(range(N),range(N))
    return qc
#
def bipartiteWalk(N,n,qreg,qcoin):
    qreg = QuantumRegister(N)
    qcoin = QuantumRegister(n)
    qc = QuantumCircuit(qreg,qcoin,name='BipartiteGraph')
    qc.x(qreg[N-1])
    qc.swap(qreg[0:N-1],qcoin[0:n])
    return qc

filePath = 'CoinedQuantumWalk/Search/'
defaultFileName = "CoinedQWSearch_N"
markedVertex = [0]
N = 3
times = 2

coinedSearchCircuit = runWalkComplete2(markedVertex,N,times)
result = simul(coinedSearchCircuit,False)
resultFig = plot_histogram(result)
savefig(resultFig,filePath,defaultFileName)
coinedSearchCircuit.draw(output='mpl')
#plt.show()

#qreg = QuantumRegister(N,'vertices')
#qcoin = QuantumRegister(N,'coin')
#creg = ClassicalRegister(N)
#qc = QuantumCircuit(qreg,qcoin,creg)

#markedVertex2=markedListComplete(markedVertex,N)
#qcOracle = oracleComplete(markedVertex2,N,False)
#qcDif = diffusionComplete(N)
#qcQWalk = completeGraphWalk(N)
#qcQWalk.draw(output='mpl')
#qc.append(qcDif,range(2*N))
#qc.append(qcOracle,range(2*N))
#qc.append(qcQWalk,range(2*N))
#qc = transpile(qc,basis_gates=['cx','u3'])
#qc.draw(output='mpl')
#plt.show()
 
