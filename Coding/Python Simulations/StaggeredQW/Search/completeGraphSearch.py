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

def completeTess(N):
    H.append((2*ones((N,N)))/N - eye(N))
    return H

def completeTessList(N):
    H=[]
    for n in N:
        H.append((2*ones((n,n)))/n - eye(n))
    return H

def oracle(N,marked):
    R = eye(N)
    R[marked][marked] = -1
    return R

def oracleList(N,marked):
    R=[]
    for n in N:
        R.append(eye(n))
    for oracle in R:
        oracle[marked][marked] = -1
    return R

def evo(theta,H,oracle):
    U = linalg.expm(-1j*theta*H)
    Uprime = U.dot(oracle)
    return Uprime

def fin(N,evo):
    psiN = init(N)
    psiN = evo.dot(psiN)
    return psiN

def plotSearch(N,probT,tSpace,configVec):
    print(tSpace)
    for steps,walk,config,n in zip(tSpace,probT,configVec,N):
        # print(config)
        # print(steps)
        plot(walk,color=config[0],linestyle=config[1],label="N=%s"%n)
        vlines(max(steps),0,1,color=config[0],linestyle=config[1])
        legend()
        xlabel("Number of steps")
        ylabel("Probability of the marked element")

def spaceGen(N):
    stepVec = []
    tVec = []
    for n in N:
        idealSteps = floor((pi/4)*sqrt(n))
        stepVec.append(int(idealSteps))
    return stepVec

def staggeredSearchList(N,tSpace,marked,oracleList,completeTessList,theta,configVec):
    prob = []
    probT = []
    pairs = []
    stepsAux = []
    stepsAuxT = []
    for n,oracle,tess,steps in zip(N,oracleList,completeTessList,tSpace):
        evol = evo(theta,tess,oracle)
        psiN = init(n)
        for step in range(1,steps+1):
            psiN = evol.dot(psiN)
            prob += [(absolute(psiN[marked][0])**2)]
            # pairs.append(((absolute(psiN[marked][0])**2),step))
            stepsAux.append(step)
        probT.append(prob)
        stepsAuxT.append(stepsAux)
        # print("Experimental steps:%s\tTheoretical Steps:%s\n"%(max(pairs),(pi/4)*sqrt(n)))
        # for obj in pairs:
            # print(obj)
        pairs = []
        stepsAux = []
        prob = []
    # print(probT)
    plotSearch(N,probT,stepsAuxT,configVec)
    show()

def staggeredSearch(N,U,steps,marked):
    psiN = init(N)
    probs = zeros((steps+1,1))
    for t in range(1,steps+1):
        psiN = U.dot(psiN)
        probAux = ampToProb(N,psiN,marked)
        probs[t] = probAux[marked]
    return psiN,probs


N=[20,40]
marked = 0
theta = pi/2
tVec = spaceGen(N)
# print(tVec)


H=completeTessList(N)
oracle = oracleList(N,marked)

colors = ['r','b','g','k']
lines = ['-','--','-.',':']
configVec = zip(colors,lines)

staggeredSearchList(N, tVec, marked, oracle, H, theta,configVec)

# print(H)
# print(oracle)
# H = completeTess(N)
# O = oracle(N,marked)
# U = completeEvo(theta,H,O)
# W,probs = staggeredSearch(N,U,steps,marked)
# probs2 = ampToProb(N,W,marked)
# print(probs2)
# print(steps)
# print("\nIdeal steps: %s\nFinal state for %s steps:\n%s"%(idealSteps,steps,W))
# print("\nAmplitudes:\n%s"%amps)

# plot(probs)
# xlabel("Steps")
# ylabel("Probability")
# show()