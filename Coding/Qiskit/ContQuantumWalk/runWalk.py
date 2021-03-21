import sys
sys.path.append('../Tools')
from IBMTools import( 
        simul,
        savefig,
        saveMultipleHist,
        printDict,
        plotMultipleQiskit,
        plotMultipleQiskitIbm,
        plotMultipleQiskitIbmSim,
        plotMultipleQiskitIbmSim2,
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
        circ.x(qreg[0])
        circ.measure(qreg,creg)
        circ = transpile(circ)
        return circ 
    else:
        circ.x(qreg[0])
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

def multDiagUniOp(N,NCirc,gamma,adjacency,time,backend,method):
    unitaryCircList = []
    for t in time:
        U0 = unitary_ctqw(gamma,N,adjacency,[],t)
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
   

filePath = 'ContQuantumWalk/'
circFilePath = 'ContQuantumWalk/Circuits/'
circDiagFilePath = 'ContQuantumWalk/Circuits/'
circQftFilePath = 'ContQuantumWalk/Circuits/'
defaultFileName = "ContQW_"
circDefaultFileName = "circContQW_"
circQftDefaultFileName = "circQft_"
circDiagDefaultFileName = "circDiag_"
style = {'figwidth':20,'fontsize':17,'subfontsize':14,'compress':True}
styleQft = {'figwidth':15,'fontsize':17,'subfontsize':14}#, 'compress':True}
styleDiag = {'figwidth':15,'fontsize':17,'subfontsize':14 }

#Cycle example
#N = 8
#c = [0,1] + [0 for x in range(N-3)] + [1]
#A = circulant_adjacency(N,c)
#G = nx.from_numpy_matrix(np.array(A))
#nx.draw(G, with_labels = True)
#plt.show()
#print('Degree: ', nx.degree(G),'\nDensity:', nx.density(G)) 

#IBMQ.load_account()
#provider = setProvider('ibm-q-minho','academicprojects','quantalab')
##leastBusyBackend =leastBusy(10,provider)
##print("Least busy backend:",leastBusyBackend)
#melBackend = provider.get_backend('ibmq_16_melbourne')
##32QV
#bogBackend = provider.get_backend('ibmq_bogota')
#parisBackend = provider.get_backend('ibmq_paris')
#manhatBackend = provider.get_backend('ibmq_manhattan')
#torontoBackend = provider.get_backend('ibmq_toronto')
#casablancaBackend = provider.get_backend('ibmq_casablanca')
##Chosen
#backend = casablancaBackend 
#simulator = provider.get_backend('ibmq_qasm_simulator')
#method = 'noise_adaptive'
method = 'trivial'
simulator = Aer.get_backend('qasm_simulator')
backend = simulator
backend2 = 'ibmq_casablanca' 
#Cont operator.##
N = 8 
NCirc = 3
gamma =  1/(2*np.sqrt(2))
t = 3
time = [0,1,2,3]
c = [0,1] + [0 for x in range(N-3)] + [1]
qft = dft(N, scale = 'sqrtn')
iqft = inv(qft)

A = iqft@circulant_adjacency(N,c)@qft
diagA = np.diag(A)
U0 = unitary_ctqw(gamma, N, A, [],t)
diagU0 = np.diag(U0).tolist()
U = iqft@U0@qft

shots = 3000
UCirc = diagUniOp(NCirc,diagU0,backend,method)
#continuousCirc = contCirc(NCirc,UCirc,backend)

unitaryCircList = multDiagUniOp(N,NCirc,gamma,A,time,backend,method)
multipleCircs = multContCirc(NCirc,unitaryCircList,time,backend,method)
#job = getJob('60509f13b24affd375a45b44',provider,backend) 
#job=backend.retrieve_job('60513dcc7240bdc91ac83cf4')
#job = execute(multipleCircs,backend=backend,shots=shots)
#job_monitor(job)
#result = job.result()
counts = [{'000': 180, '001': 2176, '010': 95, '011': 155, '100': 65, '101': 192, '110': 80, '111': 57}, {'000': 139, '001': 1670, '010': 127, '011': 198, '100': 63, '101': 386, '110': 333, '111': 84}, {'000': 100, '001': 740, '010': 257, '011': 240, '100': 131, '101': 783, '110': 595, '111': 154}, {'000': 78, '001': 184, '010': 486, '011': 338, '100': 291, '101': 700, '110': 715, '111': 208}]
#counts = result.get_counts()
#print(counts)
#plotMultipleQiskitIbm(NCirc,counts,time,shots,True)
#plt.bar(correctedResult.keys(),correctedResult.values())
#plot_histogram(correctedResult)
#g = plt.figure(1)
#g.show()

#distFig = plotMultipleQiskitIbmSim2(NCirc,multipleCircs,counts,time,shots,True,backend2)
#saveContWalkFig2(NCirc,time,distFig,filePath,defaultFileName)

#gateCountFileName = 'gateCount_N3'
#timeGateCount = np.arange(0,99,1).tolist()
#unitaryCircListGateCount = multDiagUniOp(N,NCirc,gamma,A,timeGateCount,backend,method)
#multipleCircsGateCount = multContCirc(NCirc,unitaryCircListGateCount,timeGateCount,backend,method)
#gateCountList = countTimeGates(multipleCircsGateCount)
#plotCountTimeGates(gateCountList)
#plt.savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Qiskit/'+filePath+gateCountFileName)
#   
circFigN = 3
circFigSteps = [1]
circFig = drawCirc(NCirc,UCirc,1,style)
#plt.show()
saveContWalkFig(circFigN,circFigSteps,circFig,circFilePath,circDefaultFileName)
#
circQftN = 3
circQftSteps = [1]
circQftFig = drawQftCirc(circQftN,styleQft)
saveContWalkFig(circQftN,circQftSteps,circQftFig,circQftFilePath,circQftDefaultFileName)
##
circDiagN = 3
circDiagSteps = [3]
circDiagGamma= 1/(2*np.sqrt(2))
circU0 = unitary_ctqw(circDiagGamma, N, A, [],3)
circDiagU0 = np.diag(circU0).tolist()
circDiagFig = drawDiagUni(circDiagN,circDiagU0,backend,method,styleDiag)
saveContWalkFig(circDiagN,circDiagSteps,circDiagFig,circDiagFilePath,circDiagDefaultFileName)

#initCond = '0'
#initState = init_state(N,initCond)
#psiN = final_state(U,initState)
#probvec = prob_vec(psiN,N)
#print(probvec)
#job2 = execute(multipleCircs,backend=simulator,shots=shots)
#counts2 = job2.result().get_counts()
#print(counts2)
#plotMultipleQiskitIbm(NCirc,counts2,time,shots,True)

#x = np.arange(0,N,1)
#plt.plot(x,probvec) #plot the lines
#plt.title('Caminhada Circulante')
#plt.xlabel("Graph Node")
#plt.ylabel("Probability"
#f = plt.figure(1)
#f.show()

#input()
