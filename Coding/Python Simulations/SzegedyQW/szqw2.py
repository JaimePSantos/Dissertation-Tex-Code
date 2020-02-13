from numpy import *
import matplotlib.pyplot as plt
import matplotlib
from graphs import *
from scipy import linalg
import networkx as nx
import sys
set_printoptions(threshold=sys.maxsize)

#Not working
def prob(graph,N,vertex):
    p=0
    for j in range(N):
        k = str(vertex)
        l = str(j)
        if graph.is_adjacent(k,l):
            p  =  sqrt(graph.vertex_degreeprob(k))
    return p

def graph_state(N,st):
    mat =eye(N)
    state = zeros(N)
    state = mat[:,st]
    return state

def alpha(graph,N,st):
    g = graph.adjacency_matrix()
    p = prob(graph,N,st)
    state = graph_state(N,st)
    alfa = dot(g[st],p)
    alfax = kron(state,alfa)
    return alfax

def beta(graph,N,st):
    g = graph.adjacency_matrix()
    p = prob(graph,N,st)
    state = graph_state(N,st)
    beta = dot(g[st],p)
    betay=kron(beta,state)
    return betay

def AA(graph,N):
    m = zeros(N**2)
    for i in range(0,N):
        m = m+outer(alpha(graph,N,i),alpha(graph,N,i))
    return m

def BB(graph,N):
    m = zeros(N**2)
    for i in range(0,N):
        m = m+outer(beta(graph,N,i),beta(graph,N,i))
    return m

def reflectA(graph,N):
    Ra = zeros(N**2)
    A = dot(AA(graph,N),2)
    Ra = around(A-eye(N*N))
    return Ra

def reflectB(graph,N):
    Rb = zeros(N**2)
    B = dot(BB(graph,N),2)
    Rb = around(B-eye(N**2))    
    return Rb

def operator(graph,N):
    Ra= reflectB(graph,N)
    Rb= reflectA(graph,N)
    Op = dot(Rb,Ra)
    return Op

def plotqw(N,prob):
    plot(arange(N**2),prob) 
    xlabel("Graph Node")
    ylabel("Probability")
    show()


def init_state(N,graph,st):
    psi0 = np.zeros((N**2))
    psi0[int(0)] = 1/sqrt(2)
    psi0[int((N**2)-1)] = 1/sqrt(2)

    return psi0

def final_state(Op,psi0):

    t = 100
    for i in range(t):
        psi0 = dot(Op,psi0)
   

    return psi0

def prob_vec(psiN,N):
    probs = zeros((N**2,1))
    for x in range(N**2):
        probs[x]=psiN[x]*conjugate(psiN[x])
    return probs

def szqwalk(graph,N):
    opera = operator(graph,N)
    psi0= init_state(N,graph,int(N-1))
    psiN= dot(opera,psi0)
    probvec=prob_vec(psiN,N)
    plotqw(N,probvec)


N=10
g = Graph({})
graph = g.linegraph(N)
szqwalk(graph,N)
















# print("Operator")
# print(operator(graph,N))
# print("Operator")
# print(operator(graph,N))

# e = AA(graph,N)
# f = BB(graph,N)
# print("AA")
# print(e)
# print("BB")
# print(f)
# d = beta(graph,N,0)
# print("## Alpha ## ")
# print(c)
# print("## Beta ##")
# print(d)
