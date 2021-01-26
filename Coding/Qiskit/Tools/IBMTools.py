from qiskit import IBMQ
from qiskit.tools.monitor import job_monitor
from qiskit.providers.ibmq import least_busy
from qiskit import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from qiskit.visualization import( plot_histogram,
                        plot_state_city)

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
