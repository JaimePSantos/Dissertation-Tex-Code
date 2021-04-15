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
        plotMultipleQiskitContSearch,
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

def exp_diag_qft(A,N):
    qft = dft(N, scale = 'sqrtn') #-- fourier transform
    iqft = inv(qft) #--
    D = np.diag(iqft@A@qft)
    D = np.exp(-1j*D)
    return list(D)

def diffusion_qc(expD, nq, qft_d):
    qreg = QuantumRegister(nq)
    qc = QuantumCircuit(qreg, name = 'Diagonal')
    qc.append(QFT(nq,do_swaps=False,approximation_degree = qft_d,inverse=True), range(nq))
    qc.barrier()
    qc.diagonal(expD, qreg)
    qc.barrier()
    qc.append(QFT(nq,do_swaps=False,approximation_degree = qft_d,inverse=False), range(nq))
    return qc

def oracle(N,markedList,t,r):
    O = np.zeros(N)
    for marked in markedList:
        O[marked] = 1
    O = list(np.exp(1j * O * t / r))
    return O

def oracleCirc(N,O):
    qreg = QuantumRegister(N)
    qc = QuantumCircuit(qreg,name='    Oracle    ')
    qc.diagonal(O,qreg)
    return qc

def contSearchCirc(N,NCirc,time,nTrotter,approxQFT,oracle,expD):
    qreg = QuantumRegister(NCirc)
    creg = ClassicalRegister(NCirc)
    qc = QuantumCircuit(qreg,creg)
    qcOracle = oracleCirc(NCirc,oracle)
    qcDiffusion = diffusion_qc(expD,NCirc,approxQFT)
    qc.h(qreg)
    qc.barrier()
    for n in range(nTrotter):
        qc.append(qcOracle,range(NCirc))
        qc.barrier()
        qc.append(qcDiffusion,range(NCirc))
        qc.barrier()
    qc.measure(qreg,creg)
    qc = transpile(qc,basis_gates=['cx','cp','rz','h','x'])
    return qc

def multExpD(N,A,gamma,time,r):
    expDList = []
    qft = dft(N, scale = 'sqrtn') 
    iqft = inv(qft)
    for t in time:
        B = -gamma * A * t / r
        lambdA = iqft@B@qft
        D = np.diag(lambdA)
        D = np.exp(-1j*D)
        expDList.append(list(D))
    return expDList 

def multContCircSearch(N,NCirc,expDList,time,backend,methods,approxQFT,oracle,nTrotter):
    circList = []
    qreg = QuantumRegister(NCirc)
    creg = ClassicalRegister(NCirc)
    for t,expD in zip(time,expDList):
        circ = QuantumCircuit(qreg,creg)
        circ = contSearchCirc(N,NCirc,time,nTrotter,approxQFT,oracle,expD)
        circ = transpile(circ,optimization_level=3,backend=backend, layout_method=method)
        circList.append(circ)
    return circList

def drawDiffusion(expD, nq, qft_d):
    qreg = QuantumRegister(nq)
    qc = QuantumCircuit(qreg,name='     Adj     ')
    qc.append(QFT(nq,do_swaps=False,approximation_degree = qft_d,inverse=True), range(nq))
    qc.barrier()
    qc.diagonal(expD, qreg)
    qc.barrier()
    qc.append(QFT(nq,do_swaps=False,approximation_degree = qft_d,inverse=False), range(nq))
    return qc

def drawCirc(N,NCirc,time,style,nTrotter,approxQFT,oracle,expD):
    qreg = QuantumRegister(NCirc)
    creg = ClassicalRegister(NCirc)
    qc = QuantumCircuit(qreg,creg)
    qcOracle = oracleCirc(NCirc,oracle)
    qcDiffusion = drawDiffusion(expD,NCirc,approxQFT)
    qc.h(qreg)
    qc.barrier()
    for n in range(nTrotter):
        qc.append(qcOracle,range(NCirc))
        qc.append(QFT(nq,do_swaps=False,approximation_degree = approxQFT,inverse=True,name='    QFT    '), range(nq))
        qc.append(qcDiffusion,range(NCirc))
        qc.append(QFT(nq,do_swaps=False,approximation_degree = approxQFT,inverse=False,name='    IQFT    '), range(nq))
        qc.barrier()
    qc.measure(qreg,creg)
    fig = qc.draw(output='mpl',style=style,fold=-1)
    return fig 

def drawOracleCirc(N,O,style):
    qreg = QuantumRegister(N)
    qc = QuantumCircuit(qreg,name='    Oracle    ')
    qc.diagonal(O,qreg)
    qc = transpile(qc,basis_gates=['cp','h','cx','rz'])
    fig = qc.draw(output='mpl',style=style)
    return fig

def drawAdjCirc(expD, nq, style):
    qreg = QuantumRegister(nq)
    qc = QuantumCircuit(qreg,name='     Adj     ')
    qc.diagonal(expD, qreg)
    qc = transpile(qc,basis_gates=['cp','h','cx','rz'])
    fig = qc.draw(output='mpl',style=style)
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

def runWalkResults(walkCirc,shots):
    walkResult = simul(walkCirc,False,shots)
    correctedResult = { int(k[::-1],2) : v/shots for k, v in walkResult.items()}
    return correctedResult

filePath = 'ContQuantumWalk/Search/'
circFilePath = 'ContQuantumWalk/Search/Circuits/'
defaultFileName = "ContQW_"
circDefaultFileName = "circContSearch_"
circQftDefaultFileName = "circQft_"
circAdjDefaultFileName = "circAjd_"
circOracleDefaultFileName = "circOracle_"
style = {'figwidth':20,'fontsize':17,'subfontsize':14}
styleQft = {'figwidth':15,'fontsize':17,'subfontsize':14}
styleAdj = {'figwidth':15,'fontsize':17,'subfontsize':14 }
styleOracle = {'figwidth':15,'fontsize':17,'subfontsize':14 }

backend = Aer.get_backend('qasm_simulator')
method = 'trivial'
nq = 3
N = 2 ** nq
r = 2
shots = 3000
gamma = 1 / N
t = ((np.pi/2) * np.sqrt(N))
time = [0,t/2,t,t+t/2]
time = [round(x,2) for x in time]
markedList = [1]
approxQFT = 0
cComplete = [0]+[1 for x in range(N-1)]

A = circulant_adjacency(N, cComplete)
O = oracle(N,markedList,t,r)
expD = exp_diag_qft(-gamma * A * t / r, N)
#searchCirc = contSearchCirc(N,nq,t,r,approxQFT,O,expD)
#searchResults = runWalkResults(searchCirc,shots)
#searchCirc.draw(output='mpl')
#plt.show()
#plot_histogram(searchResults)
#plt.show()

expDList = multExpD(N,A,gamma,time,r)
#print(expDList)
searchCircList = multContCircSearch(N,nq,expDList,time,backend,method,approxQFT,O,r)
#for circ in searchCircList:
#    circ.draw(output='mpl')
searchMeasFig = plotMultipleQiskitContSearch(nq,searchCircList,time,shots,True)
saveContWalkFig(nq,[r],searchMeasFig,filePath,defaultFileName)

#baseCirc = drawCirc(N,nq,t,style,r,approxQFT,O,expD)
#saveContWalkFig(nq,[r],baseCirc,circFilePath,circDefaultFileName)
#
#oracleCirc =  drawOracleCirc(nq,O,styleOracle)
#saveContWalkFig(nq,[r],oracleCirc,circFilePath,circOracleDefaultFileName)
#
#adjCirc = drawAdjCirc(expD, nq, styleAdj)
#saveContWalkFig(nq,[r],adjCirc,circFilePath,circAdjDefaultFileName)
