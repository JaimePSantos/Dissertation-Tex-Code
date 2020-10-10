from numpy import *
import matplotlib.pyplot as plt
from graphs import *
from scipy import linalg

def init_state(N): #generalizar isto ?
    psi0 = np.zeros((N,1))
    # psi0[int(N/2)-1] = 1/sqrt(2)
    # psi0[int(N/2)+1] = 1/sqrt(2)
    # psi0[0] = 1/sqrt(2)
    psi0[int(N/2)-1] = 1/sqrt(2)
    psi0[int(N/2)] = 1/sqrt(2)
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

def plotmultqw(N,prob,prob2,prob3):
    # matplotlib.rcParams.update({'font.size': 14})

    x = arange(-N/2,N/2)
    plot(x,prob,'r',label="2")
    plot(x,prob2,'g',label="4",linestyle="--")
    plot(x,prob3,'b',label="8",linestyle="-") #plot the lines
    legend(["2", "4","8"])
    xlabel("Graph Node")
    ylabel("Probability")
    show()

def plotqw(N,prob):
    matplotlib.rcParams.update({'font.size': 14})
    fig = figure()
    ax=fig.add_subplot(1,1,1) #add_subplot(nrows, ncols, index, **kwargs)
    x = np.linspace(-100,100,N)
    plot(x,prob) #plot the lines
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
    return probvec

    # plotqw(N,probvec)

N = 200
t = 100
gamma = 1/(2*np.sqrt(2)) #gamma otimo que aumenta a distancia das cristas
# gamma = 1/(2*np.sqrt(3))

lg = Graph({})
lg = lg.linegraph(N) #defined in graph.py, must have "pip install natsort" -> Future work: Remove the dependency

walkZero=ctqwalk(lg,N,t,gamma)
# walkOne= ctqwalk(lg,N,t,1/(2*np.sqrt(4)))
# walkTwo= ctqwalk(lg,N,t,1/(2*np.sqrt(8)))
# plotmultqw(N,walkZero,walkOne,walkTwo)

plotqw(N,walkZero)

