from numpy import *
import matplotlib.pyplot as plt
import matplotlib
from graphs import *
from scipy import linalg
import networkx as nx
import sys
set_printoptions(threshold=sys.maxsize)

def prob(graph,N,vertex):
    p=0
    for j in range(2*N):
        k = str(vertex)
        l = str(j)
        if graph.is_adjacent(k,l):
            p  =  p + sqrt(graph.vertex_degreeprob(k))
    return 1/p

def graph_state(graph,N,st):
    mat =eye(2*N)
    state = zeros(2*N)
    state = mat[:,st]
    return state

def alpha(graph,N,st):
    g = graph.adjacency_matrix()
    state = graph_state(graph,N,st)
    p = prob(graph,N,st)
    alfa = dot(g[st],p)
    alfax = kron(state,alfa)
    return alfax

def beta(graph,N,st):
    g = graph.adjacency_matrix()
   # print(g)
    p = prob(graph,N,st)
    state = graph_state(graph,N,st+N)
    beta = dot(g[st],p)
    betay=kron(beta,state)
    return betay

def AA(graph,N):
    m = zeros(2*N)
    for i in range(0,N):
        m = outer(alpha(graph,N,i),alpha(graph,N,i))
    return m

def BB(graph,N):
    m = zeros(2*N)
    for i in range(N,2*N):
        m = outer(beta(graph,N,i),beta(graph,N,i))
    return m

def reflectA(graph,N):
    Ra = zeros(2*N)
    A = 2*AA(graph,N)
    Ra = A - eye(2*N)
    return Ra

def reflectB(graph,N):
    Rb = zeros(2*N)
    B = 2*BB(graph,N)
    Rb = B - eye(2*N)
    return Rb

def op(Ra,Rb):
    op = dot(Ra,Rb)
    return op

def plotqw(N,prob):
    plot(arange(2*N),prob) #plot the lines
    xlabel("Graph Node")
    ylabel("Probability")
    show()


def init_state(N): #generalizar isto ?
    psi0 = np.zeros((N,1))
    psi0[int(N/2)] = 1
    return psi0

def final_state(Op,psi0):
    psiN = dot(Op,psi0)
    return psiN

def prob_vec(psiN,N):
    probs = zeros((N,1))
    for x in range(N):
        probs[x]=psiN[x]*conjugate(psiN[x]) #duvida aqui
    return probs

N=3
g = Graph({})
graph = g.bipartite_linegraph(N)
#print(alpha(graph,N,0))
#print(beta(graph,N,0))
c = alpha(graph,N,0)+alpha(graph,N,1)+alpha(graph,N,2)
d = beta(graph,N,0)+beta(graph,N,1)+beta(graph,N,2)

print(c)
print("\n")
print(d)

#print(AA(graph,N))
#print(reflectB(graph,N))
#print(BB(graph,N))
#Ra = reflectA(graph,N)
#Rb = reflectB(graph,N)
#opera=op(Ra,Rb)

#psi0 = init_state(2*N)

#psiN = final_state(opera,psi0)

#probvec=prob_vec(psiN,2*N)

#plotqw(2*N,probvec)





#print(prob(graph,N,0))

#print(alpha(graph,N,5))
#print(beta(graph,N,5))

#adj = graph.adjacency_matrix()
#state = graph_state(graph,N,1)

#print(state)




#### 
# m1 = matrix([[0,1],[0,0]])
#m2 = matrix([[0,0],[1,0]])
#mfinal = zeros((2*N,2*N))
#final = kron(m1,adj) + kron(m2,adj)

#x = array([1,0,0])

#p = sqrt(mfinal[0]/sum(mfinal[0]))

#print(mfinal[0])


#alfa0 = kron(x,p)
#print(alfa0)
#print("\n")
#print(mfinal)
####



