#CoinedQuantumWalk

from numpy import *
from matplotlib.pyplot import *
rcParams['figure.figsize'] = 11, 8
matplotlib.rcParams.update({'font.size': 15})


def coins(Matrix):
    if Matrix == "H":
        coin = array([[1/sqrt(2) , 1/sqrt(2)],[1/sqrt(2) , -1/sqrt(2)]])
    elif Matrix == "X":
        coin = array([0,1],[1,0])
    return coin

def init_state(N,Pos,Amp,CoinState):
    initstate = zeros((N,1))
    for Node in range(0,Pos.size):
        initstate[Pos[Node],0] = Amp[Node]
    g = kron(CoinState,initstate)
    return g

def init_cond(init,state0,state1):
    if init == "0":
        psi0 = array([[1],[0]])
    if init == "1":
        psi0 = array([[0],[1]])
    if init == "01":
        psi0 = array([[1/sqrt(2)],[(1*1j)/sqrt(2)]])
    return psi0

def walk_op(N,state0,state1):
    c00= outer(state0,state0)
    c01= outer(state0,state1)
    c10= outer(state1,state0)
    c11= outer(state1,state1)
    ShiftPlus = roll(eye(N),1,axis=0)
    ShiftMinus = roll(eye(N),-1,axis=0)
    Shift = kron(c00,ShiftPlus) + kron(c11,ShiftMinus) 
    return Shift

def CU_op(coin, shift,N):
    U = shift.dot(kron(coin,eye(N)))
    return U

def final_state(U,psi0,steps):
    for t in range(0,steps):
        psi0=U.dot(psi0)
    return psi0

def prob_vec(psiN,N):
    probs = zeros((N,1))
    for x in range(N):
        probs[x]=psiN[x]*conjugate(psiN[x]) + psiN[N+x]*conjugate(psiN[N+x]) #duvida aqui
    return probs

def plotqw(N,prob,init):
    x = arange(-N/2,N/2)
    plot(x,prob) 
    xlabel("Graph Node")
    ylabel("Probability")
#    savefig(r'C:\Users\Jaime\Documents\GitHub\Jaime-Santos-Dissertation\Results\Simulations\CoinedQuantumWalk\Coinedpsi0'+str(init))
    savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Simulations/CoinedQuantumWalk/Coinedpsi0'+str(init))
    clf()

def plotmultqw(N,probT,init,steps,configVec):
    x = arange(-N/2,N/2)
    plot(x,prob1,'b',label="Steps= %s"%str(steps1))
    plot(x,prob2,'g',label="Steps= %s"%str(steps2))
    plot(x,prob3,'r',label="Steps= %s"%str(steps3))
    legend()
    xlabel("Graph Node")
    ylabel("Probability")
#    savefig(r'C:\Users\Jaime\Documents\GitHub\Jaime-Santos-Dissertation\Results\Simulations\CoinedQuantumWalk\CoinedMultiplepsi001')
    savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Simulations/CoinedQuantumWalk/CoinedMultiplepsi001')
    clf()

def cqwalk(N,Steps,state0,state1,initcond):
    P = int((N+1)/2)
    Coin = coins("H")
    shift= walk_op(N,state0,state1)
    U = CU_op(Coin,shift,N)
    coinstate = initcond
    amp = array([1])
    psi0 = init_state(N,array([P]),amp,coinstate)
    psiN = final_state(U,psi0,Steps)
    probvec = prob_vec(psiN,N)
    return probvec

def multipleCoined(n,steps,state0,state1,initcond):
    probvec = []
    for step in steps:
        P = int((n+1)/2)
        Coin = coins("H")
        shift= walk_op(N,state0,state1)
        U = CU_op(Coin,shift,n)
        coinstate = initcond
        amp = array([1])
        psi0 = init_state(n,array([P]),amp,coinstate)
        psiN = final_state(U,psi0,step)
        probvec.append(prob_vec(psiN,n))
    return probvec
    
def plotmultqw2(N,probT,init,steps,configVec):
    x = arange(-N/2,N/2)
    stepsName=""
    for walk,config,step in zip(probT,configVec,steps):
        plot(x,walk,color=config[0],linestyle=config[1],label="Steps=%s"%step)
        stepsName+=str(step)
    legend()
    xlabel("Graph Node")
    ylabel("Probability")
#    savefig(r'C:\Users\Jaime\Documents\GitHub\Jaime-Santos-Dissertation\Results\Simulations\CoinedQuantumWalk\CoinedMultiplepsi001')
    savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Simulations/CoinedQuantumWalk/CoinedMultiple_psi'+str(init)+'_'+str(stepsName))
    clf()

N = 200
steps = 100
steps1 = 40
steps2 = 80
steps3 = 120

state0 = array([1,0])
state1 = array([0,1])

init = '01'
initcond = init_cond(init,state0,state1)
init1 = '1'
initcond1 = init_cond(init1,state0,state1)
init2 = '0'
initcond2 = init_cond(init2,state0,state1)

# # Single Plots
#qw = cqwalk(N,steps,state0,state1,initcond)
#plotqw(N,qw,init)
#qwx = cqwalk(N,steps,state0,state1,initcond1)
#plotqw(N,qwx,init1)
#qwy = cqwalk(N,steps,state0,state1,initcond2)
#plotqw(N,qwy,init2)
#
## # Multiple plots
#qw1=cqwalk(N,steps1,state0,state1,initcond)
#qw2=cqwalk(N,steps2,state0,state1,initcond)
#qw3=cqwalk(N,steps3,state0,state1,initcond)
#plotmultqw(N,qw1,qw2,qw3,steps1,steps2,steps3)
colors = ['r','b','g','k']
lines = ['-','-','-','-']
configVec = zip(colors,lines)

stepMult = [32,64,128]
multQWpsi0 = multipleCoined(N,stepMult,state0,state1,initcond2)
plotmultqw2(N,multQWpsi0,init2,stepMult,configVec)

configVec2 = zip(colors,lines)

stepMult2 = [32,64,128]
multQWpsi1 = multipleCoined(N,stepMult2,state0,state1,initcond1) 
plotmultqw2(N,multQWpsi1,init1,stepMult2,configVec2)

configVec3 = zip(colors,lines)

stepMult3 = [32,64,128]
multQWpsi01 = multipleCoined(N,stepMult3,state0,state1,initcond) 
plotmultqw2(N,multQWpsi01,init,stepMult3,configVec3)
