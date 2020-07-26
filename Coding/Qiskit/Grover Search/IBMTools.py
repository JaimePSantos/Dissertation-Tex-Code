from qiskit import IBMQ
from qiskit.tools.monitor import job_monitor
from qiskit.providers.ibmq import least_busy
from qiskit import *

IBMQ.load_account()

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