from numpy import *
import matplotlib.pyplot as plt
from graphs import *
from scipy import linalg
rcParams['figure.figsize'] = 11, 8
matplotlib.rcParams.update({'font.size': 15})

def init_state(N,initcond): #generalizar isto ?
    psi0 = np.zeros((N,1))
    if initcond == 'sup':
        psi0[int(N/2)-1] = 1/sqrt(2)
        psi0[int(N/2)] = 1/sqrt(2)
    if initcond== '0':
        psi0[int(N/2)] = 1
    return psi0

def final_state(Op,psi0):
    psiN = dot(Op,psi0)
    return psiN

def ct_evo(H,t,gamma):
    U = linalg.expm(-1j*gamma*H*t)
    return U

def prob_vec(psiN,N):
    probs = zeros((N,1))
    for x in range(N):
        probs[x]=psiN[x]*conjugate(psiN[x]) #duvida aqui
    return probs

def plotmultqw(N,prob1,prob2,prob3,label1,label2,label3,typeOfPlot,superposition):
    x = arange(-N/2,N/2)
    if typeOfPlot == 'Gamma':
        plot(x,prob1,'b',label="Gamma= %s"%str(label1))
        plot(x,prob2,'g',label="Gamma= %s"%str(label2))
        plot(x,prob3,'r',label="Gamma= %s"%str(label3))
    if typeOfPlot == 'Time':
        plot(x,prob1,'b',label="Time Interval= %s"%str(label1))
        plot(x,prob2,'g',label="Time Interval= %s"%str(label2))
        plot(x,prob3,'r',label="Time Interval= %s"%str(label3))
    legend()
    xlabel("Graph Node")
    ylabel("Probability")
    if typeOfPlot == 'Gamma':
        savefig(r'C:\Users\Jaime\Documents\GitHub\Jaime-Santos-Dissertation\Results\Simulations\ContQuantumWalk\ctqwMultipleGamma')
    if typeOfPlot == 'Time':
        if superposition:
            savefig(r'C:\Users\Jaime\Documents\GitHub\Jaime-Santos-Dissertation\Results\Simulations\ContQuantumWalk\ctqwMultipleTimeSuperposition')
        else:
            savefig(r'C:\Users\Jaime\Documents\GitHub\Jaime-Santos-Dissertation\Results\Simulations\ContQuantumWalk\ctqwMultipleTime')
    clf()


def plotqw(N,prob,superposition,initCond):
    x = np.linspace(-100,100,N)
    plot(x,prob) #plot the lines
    xlabel("Graph Node")
    ylabel("Probability")
    if superposition:
        savefig(r'C:\Users\Jaime\Documents\GitHub\Jaime-Santos-Dissertation\Results\Simulations\ContQuantumWalk\ctqwSingleSup')
    else:
        savefig(r'C:\Users\Jaime\Documents\GitHub\Jaime-Santos-Dissertation\Results\Simulations\ContQuantumWalk\ctqwSinglePsi'+str(initCond))
    clf()

def ctqwalk(G, N, t, gamma, initState):
    A = G.adjacency_matrix()
    D = G.degree_matrix()
    L = G.laplacian()
    psi0 = initState
    U = ct_evo(L,t,gamma)
    psiN = final_state(U,psi0)
    probvec = prob_vec(psiN,N)
    return probvec

    # plotqw(N,probvec)

N = 200
t = 100
t1 = 40
t2 = 80
t3 = 120
gamma = 1/(2*np.sqrt(2)) #gamma otimo que aumenta a distancia das cristas

gamma1 = 1/(2*np.sqrt(2)) 
gamma2 = 1/(3*np.sqrt(2)) 
gamma3 = 1/(4*np.sqrt(2)) 

lg = Graph({})
lg = lg.linegraph(N) #defined in graph.py, must have "pip install natsort" -> TODO: Remove the dependency

# # Single Plots
initCond = '0'
initState = init_state(N,initCond)
qw = ctqwalk(lg,N,t,gamma,initState)
plotqw(N,qw,False,initCond)

qwGamma = ctqwalk(lg,N,t,gamma3,initState)
plotqw(N,qwGamma,False,initCond+"LowerGamma")

initCondSup = 'sup'
initStateSup = init_state(N,initCondSup)
qwSup = ctqwalk(lg,N,t,gamma,initStateSup)
plotqw(N,qwSup,True,initCondSup)

# # Multiple Plots
qwg1 = ctqwalk(lg,N,t,gamma1,initState)
qwg2 = ctqwalk(lg,N,t,gamma2,initState)
qwg3 = ctqwalk(lg,N,t,gamma3,initState)
plotmultqw(N,qwg1,qwg2,qwg3,gamma1,gamma2,gamma3,'Gamma',False)

qwt1=ctqwalk(lg,N,t1,gamma,initState)
qwt2=ctqwalk(lg,N,t2,gamma,initState)
qwt3=ctqwalk(lg,N,t3,gamma,initState)
plotmultqw(N,qwt1,qwt2,qwt3,t1,t2,t3,'Time',False)

qwt1Sup=ctqwalk(lg,N,t1,gamma,initStateSup)
qwt2Sup=ctqwalk(lg,N,t2,gamma,initStateSup)
qwt3Sup=ctqwalk(lg,N,t3,gamma,initStateSup)
plotmultqw(N,qwt1Sup,qwt2Sup,qwt3Sup,t1,t2,t3,'Time',superposition=True)


