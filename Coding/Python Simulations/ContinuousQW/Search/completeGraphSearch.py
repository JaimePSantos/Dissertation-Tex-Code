from numpy import *
from matplotlib.pyplot import *
import matplotlib
from scipy import linalg
import networkx as nx
import sys
from numpy import kron
from numpy.core.umath import absolute

def init(N):
    psi0 = ones((N,1))/ sqrt(N)
    return psi0

def adjMatrix(N):
    adjM = ones((N,N)) - eye(N)
    return adjM

def adjMatrixList(N):
    adjM=[]
    for n in N:
        adjM.append(ones((n,n)) - eye(n))
    return adjM

def gammaList(N):
    gamma = []
    for n in N:
        gamma.append(1/n)
    return gamma

def hamiltonean(N,adjM,marked,gamma):
    H = -(gamma*adjM)
    H[marked][marked] = -1
    return H

def hamiltoneanList(N,adjM,marked,gammaList):
    H = []
    for (adjMatrix,gamma) in zip(adjM,gammaList):
        H.append(-(gamma*adjMatrix))
    for ham in H:
            ham[marked][marked] = -1
    return H

def evo(H,t):
    U = linalg.expm(-1j*H*t)
    return U

def fin(N,evo):
    psiN = init(N)
    psiN = evo.dot(psiN)
    return psiN


def ampToProb(N,psiN,marked):
    prob = zeros((N,1))
    probMarked = zeros((N,1))
    for x in range(N):
        prob[x] += (absolute(psiN[x])**2)
        probMarked[x] += (absolute(psiN[marked])**2)
    return prob,probMarked

def spaceGen(N,numOfSamples):
    stepVec = []
    tVec = []
    for n in N:
        stepVec.append((pi/2) * sqrt(n))
    for step in stepVec:
        tVec.append(linspace(0,step,num=numOfSamples))
    return tVec

def plotSearch(N,probT,tSpace,configVec):
    for T,walk,config,n in zip(tSpace,probT,configVec,N):
        print(config)
        plot(T,walk,color=config[0],linestyle=config[1],label="N=%s"%n)
        vlines(max(T),0,1,color=config[0],linestyle=config[2])
        legend()
        xlabel("Number of steps")
        ylabel("Probability of the marked element")

def runSearch(N,marked,tSpace,configVec,hamList):
    prob = []
    probT = []
    for (n,T,ham) in zip(N,tSpace,hamList):
        for t in T:
            evol = evo(ham,t)
            psiN = fin(n,evol)
            prob += [(absolute(psiN[marked][0])**2)]
        print("Sqrt(N):%s\tprob:%s\n"%(1/n,prob[0]))
        probT.append(prob)
        prob = []
    plotSearch(N,probT,tSpace,configVec)
    show()


NVec= [40,50,60]
marked = 0
gammaList = gammaList(NVec)
adjList = adjMatrixList(NVec)
hamList = hamiltoneanList(NVec,adjList,marked,gammaList)
nSamples = 50
TVec = spaceGen(NVec,nSamples)

colors = ['r','b','g','k']
lines = ['-','-','-','-']
lines2 = ['--','--','--','--']
configVec = zip(colors,lines,lines2)

runSearch(NVec,marked,TVec,configVec,hamList)
