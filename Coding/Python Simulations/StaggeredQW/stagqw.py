from numpy import *
from matplotlib.pyplot import *
from graphs import *
from scipy import linalg
import matplotlib.patches as mpatches
from fractions import Fraction

def init_state(N,init):
    psi0 = np.zeros((N,1))
    if init==0:
        psi0[int(floor(N/2))] = 1
    if init==1:     
        psi0[int(floor((N/2)+1))] = 1
    if init==10:
        psi0[int(floor(N/2))] = 1/sqrt(2)
        psi0[int(floor(N/2)+1)] = 1/sqrt(2)

    return psi0

def final_state(Op,psi0,steps):
    for i in range(steps):
        psi0 = dot(Op,psi0)
    return psi0

def ct_evo(Ha,Hb,theta):
    A = linalg.expm(1j*theta*Ha)
    B = linalg.expm(1j*theta*Hb)
    U = dot(B,A)
    # print("EVO OPERATOR")
    # print(U)
    return U

def prob_vec(psiN,N):
    probs = zeros((N,1))
    for x in range(N):
        probs[x]=psiN[x]*conjugate(psiN[x]) 
    return probs

def plotmultqw(N,prob,prob2,prob3,theta1,theta2,theta3):
    # matplotlib.rcParams.update({'font.size': 14})

    x = arange(-N/2,N/2)
    plot(x,prob,'b',label="pi/3")
    plot(x,prob2,'g',label="pi/4")
    plot(x,prob3,'r',label="pi/5") #plot the lines
    legend(["pi/3", "pi/4","pi/5"])
    xlabel("Graph Node")
    ylabel("Probability")
    show()

def plotqw(N,prob):
    #x = arange(N)
    x = linspace(N/2,-N/2,N)
    plot(x,prob,'r')
    xlabel("Graph Node")
    ylabel("Probability")
    show()

def stgqwalk(t1,t2, N, theta,steps,init):
    A = t1.adjacency_matrix()
    B = t2.adjacency_matrix()
    # print("Alfa")
    # print(A)
    # print("#######")
    # print("Beta")
    # print(B)
    # print("#######")
    
    Ha = A 
    Hb = B 
    psi0 = init_state(N,init)
    U = ct_evo(Ha,Hb,theta)
    psiN = final_state(U,psi0,steps)
    probvec = prob_vec(psiN,N)
    #plotqw(N,probvec)
    return probvec

N = 200
theta1= pi/3
theta2= pi/4
theta3= pi/5
steps = 50
i = Graph({})
j = Graph({})

t1 = i.line_tesselation(N,'alpha') 
t2 = j.line_tesselation(N,'beta')

#print(t1.adjacency_matrix())
#print()
#print(t2.adjacency_matrix())
k1=stgqwalk(t1,t2,N,theta1,steps,1)
k2=stgqwalk(t1,t2,N,theta2,steps,0)
k3=stgqwalk(t1,t2,N,theta3,steps,0)

# plotmultqw(N,k1,k2,k3,theta1,theta2,theta3)
plotqw(N,k1)