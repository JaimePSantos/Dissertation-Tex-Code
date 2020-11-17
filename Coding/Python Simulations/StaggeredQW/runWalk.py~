#StaggeredQuantumWalk

from numpy import *
from matplotlib.pyplot import *
from graphs import *
from scipy import linalg
import matplotlib.patches as mpatches
from fractions import Fraction
rcParams['figure.figsize'] = 11, 8
matplotlib.rcParams.update({'font.size': 15})

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

def plotmultqw(N,prob1,prob2,prob3,label1,label2,label3,initCond):
    # matplotlib.rcParams.update({'font.size': 14})

    x = arange(-N/2,N/2)
    plot(x,prob1,'b',label=r"$\theta = \frac{\pi}{%s}$"%label1)
    plot(x,prob2,'g',label=r"$\theta = \frac{\pi}{%s}$"%label2)
    plot(x,prob3,'r',label=r"$\theta = \frac{\pi}{%s}$"%label3) #plot the lines
    legend()
    xlabel("Graph Node")
    ylabel("Probability")
    # show()
#    savefig(r'C:\Users\Jaime\Documents\GitHub\Jaime-Santos-Dissertation\Results\Simulations\StagQuantumWalk\stagqwMultiple')
    savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Simulations/CoinedQuantumWalk/stagqwMultiple')
    clf()

def plotqw(N,prob,initcond):
    #x = arange(N)
    x = linspace(N/2,-N/2,N)
    plot(x,prob)
    xlabel("Graph Node")
    ylabel("Probability")
    # show()
#    savefig(r'C:\Users\Jaime\Documents\GitHub\Jaime-Santos-Dissertation\Results\Simulations\StagQuantumWalk\stagqwSingle'+str(initcond))
    savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Simulations/CoinedQuantumWalk/stagqwSingle'+str(initcond))
    clf()

def stgqwalk(t1,t2, N, theta,steps,init):
    A = t1.adjacency_matrix()
    B = t2.adjacency_matrix()    
    Ha = A 
    Hb = B 
    psi0 = init_state(N,init)
    U = ct_evo(Ha,Hb,theta)
    psiN = final_state(U,psi0,steps)
    probvec = prob_vec(psiN,N)
    #plotqw(N,probvec)
    return probvec

N = 200
denom1 = 3
denom2 = 4
denom3 = 5
initcond = 10

theta1= pi/denom1
theta2= pi/denom2
theta3= pi/denom3
steps = 50
i = Graph({})
j = Graph({})

t1 = i.line_tesselation(N,'alpha') 
t2 = j.line_tesselation(N,'beta')

qwSup1=stgqwalk(t1,t2,N,theta1,steps,initcond)
qwSup2=stgqwalk(t1,t2,N,theta2,steps,initcond)
qwSup3=stgqwalk(t1,t2,N,theta3,steps,initcond)
plotmultqw(N,qwSup1,qwSup2,qwSup3,denom1,denom2,denom3,str(initcond))

steps1 = 50
initcond0 = 0
qwSingle0=stgqwalk(t1,t2,N,theta1,steps1,initcond0)
plotqw(N,qwSingle0,initcond0)

initcond1 = 1
qwSingle1=stgqwalk(t1,t2,N,theta1,steps1,initcond1)
plotqw(N,qwSingle1,initcond1)
