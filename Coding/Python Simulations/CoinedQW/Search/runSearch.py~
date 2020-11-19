#CoinedQuantumWalkSearch

from numpy import *
from matplotlib.pyplot import *
import matplotlib
from scipy import linalg
import sys
from numpy import kron
from numpy.core.umath import absolute
matplotlib.rcParams.update({'font.size': 14})
rcParams['figure.figsize'] = 11, 8

def ket2pos(pos,dim):
    position = 0
    for i in range(len(dim)):
        position += pos[i]*prod(dim[i+1::])
    return int(position)

def init(N):
    psi0 = ones((N**2,1))/ N
    return psi0

def flipFlop(N):
    S = zeros((N**2, N**2))
    dim = [N,N]
    for v1 in range(N):
        for v2 in range(N):
            row = ket2pos([v1,v2],dim)
            column = ket2pos([v2,v1],dim)
            S[row][column] = 1
    return S

def flipFlopList(N):
    shiftList = []
    for n in N:
        S = zeros((n**2, n**2))
        dim = [n,n]
        for v1 in range(n):
            for v2 in range(n):
                row = ket2pos([v1,v2],dim)
                column = ket2pos([v2,v1],dim)
                S[row][column] = 1
        shiftList.append(S)
    return shiftList


def oracle(N,marked):
    R = eye(N**2)
    for v in range(N):
        vec = ket2pos([marked,v],[N,N])
        R[vec][vec] = R[vec][vec] - 2
    return R

def oracleList(N,marked):
    RList = []
    for n in N:
        R = eye(n**2)
        dim = [n,n]
        for v in range(n):
            vec = ket2pos([marked,v],dim)
            R[vec][vec] = R[vec][vec] - 2
        RList.append(R)
    return RList

def coin(N):
    s = 2*ones((N,N)) / N
    coin = s - eye(N)
    return coin

def coinList(N):
    coinList = []
    for n in N:
        s = 2*ones((n,n)) / n
        coin = s - eye(n)
        coinList.append(coin)
        coin=[]
        s=[]
    return coinList

def evo(N,marked):
    op = flipFlop(N).dot(kron(eye(N),coin(N))).dot(oracle(N,marked))
    return op

def evo2(N,marked,flipFlop,coin,oracle):
    coinOperator = kron(eye(N),coin)
    op = flipFlop.dot(coinOperator).dot(oracle)
    return op

def ampToProb(N,psiN):
    prob = zeros((N,1))
    for x in range(N):
        for c in range(N):
            prob[x]+= (absolute(psiN[ket2pos([x,c],[N,N])]))**2
    return prob

def spaceGen(N):
    stepVec = []
    tVec = []
    for n in N:
        idealSteps = floor((pi/2)*sqrt(n))
        stepVec.append(int(idealSteps))
    return stepVec

def plotSearch(N,probT,steps,configVec):
    plotName = ""
    stepAux = []
    stepList = []
    for step in steps:
        for i in range(1,step+1):
            stepAux.append(i)
        stepList.append(stepAux)
        stepAux = []
    for steps,walk,config,n in zip(stepList,probT,configVec,N):
        # print(walk)
        plot(walk,color=config[0],linestyle=config[1],label="N=%s"%n)
        vlines(max(steps),0,walk[-1],color=config[0],linestyle=config[2])
        legend()
        xlabel("Number of steps")
        ylabel("Probability of the marked element")
    for n in N:
        plotName+=str(n)
    savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Simulations/CoinedQuantumWalk/Search/CoinedSearch'+plotName)
    clf()

def coinedSearchList(N,tSpace,marked,shiftList,oracleList,coinList,configVec):
    prob=[]
    probAux = []
    probT=[]
    stepsAux = []
    stepsAuxT = []
    for n,shift,oracle,steps,coin in zip(N,shiftList,oracleList,tSpace,coinList):
        evol = evo2(n,marked,shift,coin,oracle)
        psiN = init(n)
        prob +=  [ampToProb(n,psiN)[marked][0]]
        # probs[0] = ampToProb(n,psiN)[marked]
        for step in range(1,steps+1):
            psiN = evol.dot(psiN)
            probAux = ampToProb(n,psiN)
            # print(probAux)
            prob += [probAux[marked][0]]
            # print(prob)
            # prob += [(absolute(psiN[marked][0])**2)]
            stepsAux.append(step)
            # print("prob dos marcados apos %s steps\n%s"%(step,prob))

        # print(prob)
        probT.append(prob)
        stepsAuxT.append(stepsAux)
        prob = []
        stepsAux = []
        probAux = []
    # print(probT)
    # for walk in probT:
        # plot(walk)
    return probT

def fin(N,marked,steps):
    psiN = init(N)
    ev = evo(N,marked)
    probs = zeros((steps+1,1))
    probs[0] = ampToProb(N,psiN)[marked]
    for i in range(1,steps+1):
        psiN = ev.dot(psiN)
        probAux = ampToProb(N,psiN)
        probs[i] = probAux[marked]
    return psiN,probs


N = [16,32,64]
marked = 0

shiftList = flipFlopList(N)
# print(shiftList)
oracleList = oracleList(N,marked)
# print(oracleList)
coinList = coinList(N)
# print(coinList)
tVec = spaceGen(N)
# print(tVec)

colors = ['r','b','g','k']
lines = ['-','-','-','-']
lines2 = ['--','--','--','--']
configVec = zip(colors,lines,lines2)

coinedSearch=coinedSearchList(N,tVec,marked,shiftList,oracleList,coinList,configVec)
plotSearch(N,coinedSearch,tVec,configVec)
