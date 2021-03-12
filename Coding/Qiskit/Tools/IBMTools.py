from qiskit import IBMQ
from qiskit.tools.monitor import job_monitor
from qiskit.providers.ibmq import least_busy
from qiskit import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from qiskit.visualization import( plot_histogram,
                        plot_state_city)
import numpy as np
#IBMQ.load_account()

def run(circuit, backend, **kwargs):
    if type(backend) is str:
        backend = Aer.get_backend(backend)
    return execute(circuit, backend, **kwargs)

def textResults(results,collisions):
    for key in results:
        if(results[key]>collisions):
            text= str(key)+ '->'+  str(results[key])
    return text

def setProvider(hub,group,project):
    provider = IBMQ.get_provider(hub=hub, group=group, project=project)
    return provider

def leastBusy(minQubits,provider):
    large_enough_devices = provider.backends(filters=lambda x: x.configuration().n_qubits > minQubits  and not x.configuration().simulator)
    leastBusybackend = least_busy(large_enough_devices)
    return leastBusybackend

def getJobCounts(result,backend):
    jobID = result.job_id()
    job = backend.retrieve_job(jobID)
    resultCount = job.result().get_counts()
    return resultCount

def listBackends(provider):
    for backend in provider.backends():
        print( backend.name())

def getJob(jobID,provider,backend):
    job = backend.retrieve_job(jobID)
    resultCount = job.result().get_counts()
    return resultCount

def printBestSeed(qc,basisGatesD,deviceBackend,startSeed,endSeed):
    dict = {}
    dict2 = {}
    for i in range(startSeed,endSeed):
        qCirc = transpile(qc,basis_gates=basisGatesD,backend=deviceBackend,optimization_level=3,layout_method='noise_adaptive',seed_transpiler=i)
        dict[i] = qCirc.count_ops()['cx']
        dict2[i] = qCirc.depth()
    print(min(dict.items(), key=lambda x: x[1])) 
    print(min(dict2.items(), key=lambda x: x[1]))

def simul(qc,stateVec,shots):
    if stateVec:
        backend = Aer.get_backend('statevector_simulator')
        result = execute(qc,backend,shots=shots).result().get_statevector(qc,decimals=3)
    else:
        backend = Aer.get_backend('qasm_simulator')
        result = execute(qc,backend,shots=3000).result().get_counts()
    return result

def savefig(fig,filePath,fileName):
    plt.savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Qiskit/'+filePath+fileName)
    plt.clf()

def saveMultipleHist(N,steps,circListList,filePath,defaultFileName):
    "Saves len(steps) histogram plots of N qubit quantum walks as .png files. Default path is /Coding/Qiskit. Default file name should be \"name_N\" "
    fileName = ""
    fileNameAux = ""
    for n,circList in zip(N,circListList):
        fileName += defaultFileName+str(n)
        for circ,step in zip(circList,steps):
            fileNameAux = "_S"+str(step)  
            result = simul(circ,False)
            resultFig = plot_histogram(result)
            fileName += fileNameAux
            fileNameAux = ""
            savefig(resultFig,filePath,fileName)
            print("%s was saved to %s"%(fileName,filePath))
            fileName = defaultFileName+str(n)
        fileName = ""

def printDict(dictionary):
    for i,k in zip(dictionary.keys(),dictionary.values()):
        print("%s: %s"%(i,k))

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

def decResultDict(n):
    "Retuns a dictionary composed of a range of N keys converted to binary."
    baseDict = {}
    for decNumber in range(2**n):
        dec = decNumber 
        baseDict[str(dec)] = 0
    return baseDict

def multDecResultDict(N,steps):
    "Returns multiple binary dictionaries."
    baseResultDictList = []
    for n in N:
        for step in steps:
            baseDict = decResultDict(n)
            baseResultDictList.append(baseDict)
    return baseResultDictList

def binResultDict(n):
    "Retuns a dictionary composed of a range of N keys converted to binary."
    baseDict = {}
    for decNumber in range(2**n):
        decToBin = bin(decNumber)[2:].zfill(ceil(log(2**n,2)))
        baseDict[str(decToBin)] = 0
    return baseDict

def multBinResultDict(N,steps):
    "Returns multiple binary dictionaries."
    baseResultDictList = []
    for n in N:
        for step in steps:
            baseDict = binResultDict(n)
            baseResultDictList.append(baseDict)
    return baseResultDictList

def multNormalizedResultDict(baseDictList,qiskitDictList):
    "Returns the result of merging qiskit produced dictionaries with dictionaries produced from multBinResultDict for graph formatting reasons."
    normalizedResultDictList = []
    for baseDict,qiskitDict in zip(baseDictList,qiskitDictList):
        normalizedResultDict = {**baseDict,**qiskitDict}
        normalizedResultDictList.append(normalizedResultDict)
    return normalizedResultDictList

def multResultsSim(multipleCircs,shots,Decimal):
    "Returns the dictionary produced by QASM simulator with the MSB changed to convention, and values (previously frequencies) converted to probabilities."
    resultList = []
    result = {}
    correctedResult = {}
    for circList in multipleCircs:
        for circ in circList:
            result = simul(circ,False,shots)
            if Decimal:
                correctedResult = { str(int(k[::-1],2)) : v/shots for k, v in result.items()}
            else:
                correctedResult = { k[::-1] : v/shots for k, v in result.items()}
            resultList.append(correctedResult)
            result = {}
    return resultList

#TODO: Delegar formatacao para uma funcao propria.
#TODO: Os labels dos eixos nao estao perfeitamente centrados. O do y fica no ultimo subplot, por alguma razao.
def multSubPlot(resultList,steps):
    "Produces a matplotlib figure composed of several subplots for different numbers of graph nodes and circuit iterations."
    nrows = len(resultList) 
    ncols = 1
    index = 1
    fig = plt.figure()
    axList = []
    auxList = []
    for resultAux,step in zip(resultList,steps):
        axList.append(fig.add_subplot(nrows,ncols,index))
        axList[-1].bar(resultAux.keys(),resultAux.values(),width=0.4,label = "Steps=%s"%step)
        axList[-1].legend()
        index+=1
    for ax in axList:
        axList[-1].get_shared_y_axes().join(axList[-1],ax)
    for ax in axList[:-1]:
        ax.set_xticklabels([])
    axList[-1].set_xticklabels(resultList[-1].keys(),rotation=45)
    plt.xlabel("Graph Node")
    plt.ylabel("Probability")
    fig.tight_layout(pad=1.0)
    return axList 

def plotMultipleQiskit(N,multipleCircs,steps,shots,Decimal):
    "Brings every dictionar and plot building functions together to either show or save the matplotlib figure."
    qiskitResultList = multResultsSim(multipleCircs,shots,Decimal)
    if Decimal:
        baseDictList = multDecResultDict(N,steps)
    else:
        baseDictList = multBinResultDict(N,steps)
    normalizedResultDictList = multNormalizedResultDict(baseDictList,qiskitResultList)
    fig = multSubPlot(normalizedResultDictList,steps)
    return fig
