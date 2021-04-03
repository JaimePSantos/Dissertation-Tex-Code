import sys
sys.path.append('../../Tools')
from IBMTools import( 
        simul,
        savefig,
        saveMultipleHist,
        printDict,
        plotMultipleQiskit,
        plotMultipleQiskitIbm,
        plotMultipleQiskitIbmSim,
        plotMultipleQiskitIbmSim2,
        plotMultipleQiskitGrover,
        plotMultipleQiskitGrover2,
        multResultsSim,
        setProvider,
        leastBusy,
        listBackends,
        getJob)
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
import networkx as nx
mpl.rcParams['figure.figsize'] = 11,8
mpl.rcParams.update({'font.size' : 15})

def circulant_adjacency(n,v): #--- it computes an adjacency matrix for the circulant graph
    iv = list(range(0,n))
    av = list(range(0,n-1))
    C = np.zeros([n,n])
    for z in range(n):
        C[z,0] = v[iv[z]]
    for x in range(1,n):
        av = iv[0:-1]
        iv[0] = iv[-1]
        iv[1::] = av
        for y in range(0,n):
            C[y,x] = v[iv[y]]
    return C

def unitary_ctqw(gamma, N, A, marked, t): #---
    Oracle = np.zeros([N,N])
    for x in marked:
        Oracle[x,x] = 1
    U = expm(1j*(-gamma*A - Oracle)*t)
    return U

def trotter(gamma, N, A, marked, t, n_trotter):
    O = np.zeros([N,N])
    for x in marked:
        O[x,x] = 1
    U = matrix_power(expm(1j*(-gamma*A)*t/n_trotter)@expm(1j*(- O)*t/n_trotter), n_trotter)
    return U

def init_state(N,initcond): #generalizar isto ?
    psi0 = np.zeros((N,1))
    if initcond == 'sup':
        psi0[int(N/2)-1] = 1/sqrt(2)
        psi0[int(N/2)] = 1/sqrt(2)
    if initcond== '0':
        psi0[0] = 1
    return psi0

def final_state(Op,psi0):
    psiN = np.dot(Op,psi0)
    return psiN

def prob_vec(psiN,N):
    probs = np.zeros((N,1))
    for x in range(N):
        probs[x]=psiN[x]*np.conjugate(psiN[x]) 
    return probs

def diagUniOp(N,diagU0,backend,method):
    qreg = QuantumRegister(N)
    creg = ClassicalRegister(N)
    circ = QuantumCircuit(qreg,name='    UniOp    ')
    circ.diagonal(diagU0,qreg) 
    circ = transpile(circ)#,optimization_level=3)#,backend=backend,layout_method=method) 
    return circ

def contCirc(N,diagUniOp,backend,method,t):
    qreg = QuantumRegister(N)
    creg = ClassicalRegister(N)
    circ = QuantumCircuit(qreg,creg)
    if t == 0: 
        #circ.h(qreg)
        circ.x(qreg[2])
        circ.measure(qreg,creg)
        circ = transpile(circ)
        return circ 
    else:
        #circ.h(qreg)
        circ.x(qreg[2])
        circ.append(QFT(N,do_swaps=False,approximation_degree=0,inverse=False,name='    QFT    '),range(N))
        circ.append(diagUniOp,range(N))
        circ.append(QFT(N,do_swaps=False,approximation_degree=0,inverse=True,name='    IQFT'    ),range(N))
        circ.measure(qreg,creg)
        circ=transpile(circ)
    return circ

def contCirc2(N,diagUniOp,backend,method,t):
    qreg = QuantumRegister(N)
    creg = ClassicalRegister(N)
    circ = QuantumCircuit(qreg,creg)
    if t == 0: 
        circ.x(qreg[0])
        circ.measure(qreg,creg)
        circ = transpile(circ)
        return circ 
    else:
        circ.append(diagUniOp,range(N))
        circ=transpile(circ)
    return circ

def multDiagUniOp(N,NCirc,gamma,adjacency,time,backend,method,marked):
    unitaryCircList = []
    for t in time:
        U0 = unitary_ctqw(gamma,N,adjacency,marked,t)
        diagU0 = np.diag(U0).tolist()
        diagCirc = diagUniOp(NCirc,diagU0,backend,method)
        unitaryCircList.append(diagCirc)
    return unitaryCircList

def multDiagUniOpTrotter(N,NCirc,gamma,adjacency,time,backend,method,marked,nTrotter):
    unitaryCircList = []
    for t in time:
        U0 = trotter(gamma, N, adjacency, marked, t, nTrotter)
        diagU0 = np.diag(U0).tolist()
        diagCirc = diagUniOp(NCirc,diagU0,backend,method)
        unitaryCircList.append(diagCirc)
    return unitaryCircList

def multContCirc(N,unitaryList,time,backend,methods):
    circList = []
    circListAux = []
    qreg = QuantumRegister(N)
    qsub = QuantumRegister(1)
    creg = ClassicalRegister(N)
    for t,diagU0 in zip(time,unitaryList):
        circ = QuantumCircuit(qreg,creg)
        circ =  contCirc(N,diagU0,backend,method,t)
        circ = transpile(circ,optimization_level=3,backend=backend, layout_method=method)
        circList.append(circ)
    return circList

def multContCirc2(N,unitaryList,time,backend,methods):
    circList = []
    circListAux = []
    qreg = QuantumRegister(N)
    qsub = QuantumRegister(1)
    creg = ClassicalRegister(N)
    for t,diagU0 in zip(time,unitaryList):
        circ = QuantumCircuit(qreg,creg)
        circ =  contCirc2(N,diagU0,backend,method,t)
        circ = transpile(circ,optimization_level=3,backend=backend, layout_method=method)
        circList.append(circ)
    return circList

def drawCirc(N,diagU,time,style):
    qreg = QuantumRegister(N)
    creg = ClassicalRegister(N)
    circ = QuantumCircuit(qreg,creg)
    circ.x(qreg[0])
    circ.barrier()
    circ.append(QFT(N,do_swaps=False,approximation_degree=0,inverse=False,name='    QFT    '),range(N))
    circ.append(diagU,range(N))
    circ.append(QFT(N,do_swaps=False,approximation_degree=0,inverse=True,name='    IQFT    '),range(N))
    circ.barrier()
    circ.measure(qreg,creg)
    fig = circ.draw(output='mpl',style=style)
    return fig 

def drawCirc2(N,diagU,time,style):
    qreg = QuantumRegister(N)
    creg = ClassicalRegister(N)
    circ = QuantumCircuit(qreg,creg)
    circ.x(qreg[0])
    circ.barrier()
    circ.append(QFT(N,do_swaps=False,approximation_degree=0,inverse=False,name='    QFT    '),range(N))
    circ.append(diagU,range(N))
    circ.append(QFT(N,do_swaps=False,approximation_degree=0,inverse=True,name='    IQFT    '),range(N))
    circ.barrier()
    circ.measure(qreg,creg)
    #fig = circ.draw(output='mpl',style=style)
    return circ 

def drawQftCirc(N,style):
    qreg = QuantumRegister(N)
    creg = ClassicalRegister(N)
    circ = QuantumCircuit(qreg,creg)
    circ = QFT(N,do_swaps=False,approximation_degree=0,inverse=False,name='    QFT    ')
    circ = transpile(circ, basis_gates=['cp','h','cx','rz'])
    fig = circ.draw(output='mpl',style=style)
    return fig 

def drawDiagUni(N,diagU0,backend,method,style):
    qreg = QuantumRegister(N)
    creg = ClassicalRegister(N)
    circ = QuantumCircuit(qreg,name='    UniOp    ')
    print(diagU0)
    circ.diagonal(diagU0,qreg) 
    #circ = transpile(circ,backend=backend,optimization_level=1)#,backend=backend,layout_method=method)#
    circ = transpile(circ,basis_gates=['cp','h','cx','rz'])
    #circ.decompose()
    fig = circ.draw(output='mpl',style=style)
    return fig

def saveContWalkFig(N,steps,fig, filePath, defaultFileName):
    specificFileName = ""
    i=0
    specificFileName+= "N%s_S"%N
    for step in steps:
        specificFileName+="%s"%step
    savefig(fig, filePath,defaultFileName+specificFileName)
    plt.clf()
    return specificFileName

def saveContWalkFig2(N,steps,fig, filePath, defaultFileName):
    mpl.rcParams.update(mpl.rcParamsDefault)
    mpl.rcParams['figure.figsize'] = 11,8
    mpl.rcParams.update({'font.size' : 15})
    specificFileName = ""
    i=0
    specificFileName+= "N%s_S"%N
    for step in steps:
        specificFileName+="%s"%step
    savefig(fig, filePath,defaultFileName+specificFileName)
    plt.clf()
    return specificFileName

def countTimeGates(circList):
    gateCountList = []
    for circ in circList:
        gateCount = circ.count_ops()
        gateCountList.append(gateCount)
        gateCount = 0
    gateCountList.pop(0)
    return gateCountList

def plotCountTimeGates(gateCountList):
    cxCount = []
    rzCount = []
    xCount = []
    cpCount = []
    for dictionary in gateCountList:
        for gate,count in dictionary.items():
            if gate == 'cx':
                cxCount.append(count)
            if gate == 'rz':
                rzCount.append(count)
            if gate == 'x':
                xCount.append(count)
            if gate == 'cp':
                cpCount.append(count)

    timeGateCount2 = np.arange(0,98,1).tolist()
    plt.plot(timeGateCount2,cxCount,label='CNot number',color='r')
    plt.plot(timeGateCount2,rzCount,label='RZ number',color='b')
    plt.yticks(np.arange(0,11,1))
    plt.xlabel("Time")
    plt.ylabel("Number of gates")
    plt.legend()
   
def trotter(gamma, N, A, marked, t, n_trotter):
    O = np.zeros([N,N])
    for x in marked:
        O[x,x] = 1
    U = matrix_power(expm(1j*(-gamma*A)*t/n_trotter)@expm(1j*(- O)*t/n_trotter), n_trotter)
    return U

filePath = 'ContQuantumWalk/'
circFilePath = 'ContQuantumWalk/Circuits/'
circDiagFilePath = 'ContQuantumWalk/Circuits/'
circQftFilePath = 'ContQuantumWalk/Circuits/'
defaultFileName = "ContQW_"
circDefaultFileName = "circContQW_"
circQftDefaultFileName = "circQft_"
circDiagDefaultFileName = "circDiag_"
style = {'figwidth':20,'fontsize':17,'subfontsize':14}#,'compress':True}
styleQft = {'figwidth':15,'fontsize':17,'subfontsize':14}#, 'compress':True}
styleDiag = {'figwidth':15,'fontsize':17,'subfontsize':14 }

backend = Aer.get_backend('qasm_simulator')
method = 'trivial'
shots = 3000
N = 8
NCirc = 3
optimalTime = (np.pi/4)*np.sqrt(N)
time = optimalTime
walkTime = 1
marked = [0]
gamma = 1 / N
walkGamma = 1 / 2*(np.sqrt(2))
cComplete = [0]+[1 for x in range(N-1)]
cCycle = [0,1] + [0 for x in range(N-3)] + [1]
qft = dft(N,scale = 'sqrtn')
iqft = inv(qft)

A = circulant_adjacency(N,cCycle)
lambdA =  iqft@A@qft

walkU0 = unitary_ctqw(walkGamma, N, lambdA, [], walkTime)
walkU = np.diag(walkU0).tolist()

walkUQiskit = diagUniOp(NCirc,walkU,backend,method)
walkCirc = contCirc(NCirc,walkUQiskit,backend,method,walkTime)
walkResult = simul(walkCirc,False,shots)
print(walkResult)
correctedResult = { int(k[::-1],2) : v for k, v in walkResult.items()}
print(correctedResult)
plot_histogram(correctedResult)
plt.show()
