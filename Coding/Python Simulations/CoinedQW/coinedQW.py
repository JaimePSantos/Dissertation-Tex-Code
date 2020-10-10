from numpy import *
from matplotlib.pyplot import *


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
    print("Shift Plus:\n",ShiftPlus)
    ShiftMinus = roll(eye(N),-1,axis=0)
    print("Shift Minus:\n",ShiftMinus)
    
    Shift = kron(c00,ShiftPlus) + kron(c11,ShiftMinus)
 
    print("Shift Matrix:\n",Shift)
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


def plotqw(N,prob):
    matplotlib.rcParams.update({'font.size': 14})
    fig = figure(figsize=(11,8))
    ax=fig.add_subplot(1,1,1) #add_subplot(nrows, ncols, index, **kwargs)
    x = range(N)
    # fig=figure(figsize=(50, 50))
    plot(x,prob) #plot the lines
    xlabel("Graph Node")
    ylabel("Probability")

    show()

def plotmultqw(N,prob,prob2,prob3,steps1,steps2,steps3):
    matplotlib.rcParams.update({'font.size': 18})
    fig = figure(figsize=(13,8))
    x = arange(-N/2,N/2)
    plot(x,prob,'b',label=str(steps1))
    plot(x,prob2,'g',label=str(steps2))
    plot(x,prob3,'r',label=str(steps3)) #plot the lines
    legend([str(steps1), str(steps2),str(steps3)])
    xlabel("Graph Node")
    ylabel("Probability")
    show()

def cqwalk(N,Steps,state0,state1):
    P = int((N+1)/2)
    
    Coin = coins("H")

    shift= walk_op(N,state0,state1)

    U = CU_op(Coin,shift,N)

    coinstate = init_cond("0",state0,state1)

    amp = array([1])
    psi0 = init_state(N,array([P]),amp,coinstate)
    psiN = final_state(U,psi0,Steps)
    probvec = prob_vec(psiN,N)
    
    return probvec
    # plotqw(N,probvec)

    
N = 200
steps1 = 40
steps2 = 80
steps3 = 120

##circuit
state0 = array([1,0])
state1 = array([0,1])
##circuit

qw1=cqwalk(N,steps1,state0,state1)
qw2=cqwalk(N,steps2,state0,state1)
qw3=cqwalk(N,steps3,state0,state1)

plotmultqw(N,qw1,qw2,qw3,steps1,steps2,steps3)