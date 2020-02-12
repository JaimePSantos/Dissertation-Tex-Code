from numpy import *
import matplotlib.pyplot as plt
from graphs import *
from scipy import linalg

def init_state(N):
    psi0 = np.zeros((N,1))
    psi0[int(0)] = 1/sqrt(2)
    psi0[int(N-1)] = 1/sqrt(2)
    return psi0

def final_state(Op,psi0,steps):
    for i in range(steps):
        psi0 = dot(Op,psi0)
    return psi0

def ct_evo(Ha,Hb,theta):
    A = linalg.expm(1j*theta*Ha)
    B = linalg.expm(1j*theta*Hb)
    U = dot(B,A)
    return U

def prob_vec(psiN,N):
    probs = zeros((N,1))
    for x in range(N):
        probs[x]=psiN[x]*conjugate(psiN[x]) 
    return probs

def plotqw(N,prob):
    plot(arange(N),prob) 
    xlabel("Graph Node")
    ylabel("Probability")
    show()

def stgqwalk(t1,t2, N, theta,steps):
    A = t1.adjacency_matrix()
    B = t2.adjacency_matrix()
    Ha = A 
    Hb = B 
    psi0 = init_state(N)
    U = ct_evo(Ha,Hb,theta)
    psiN = final_state(U,psi0,steps)
    probvec = prob_vec(psiN,N)
    plotqw(N,probvec)

N = 200
theta= pi/3
steps = 100
i = Graph({})
j = Graph({})

t1 = i.line_tesselation(N,'alpha') 
t2 = j.line_tesselation(N,'beta')

stgqwalk(t1,t2,N,theta,steps)