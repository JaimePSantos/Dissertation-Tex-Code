from numpy import *
import matplotlib.pyplot as plt
from graphs import *
from scipy import linalg

def init_state(N):
    psi0 = np.zeros((N,1))
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
        probs[x]=psiN[x]*conjugate(psiN[x]) 
    return probs

def plotqw(N,prob):
    matplotlib.rcParams.update({'font.size': 14})
    fig = figure()
    ax=fig.add_subplot(1,1,1) 
    x = np.linspace(-100,100,N)
    plot(x,prob) 
    xlabel("Graph Node")
    ylabel("Probability")
    show()

def ctqwalk(G, N, t, gamma):
    A = G.adjacency_matrix()
    D = G.degree_matrix()
    L = G.laplacian()
    psi0 = init_state(N)

    U = ct_evo(L,t,gamma)

    psiN = final_state(U,psi0)

    probvec = prob_vec(psiN,N)

    plotqw(N,probvec)

N = 200
t = 100
gamma = 1/(2*np.sqrt(2)) #gamma otimo que aumenta a distancia das cristas

lg = Graph({})
lg = lg.linegraph(N) #defined in graph.py, must have "pip install natsort" -> Future work: Remove the dependency

ctqwalk(lg,N,t,gamma)

