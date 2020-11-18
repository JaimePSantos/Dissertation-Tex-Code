#StaggeredQuantumWalkSearch

from numpy import *
from matplotlib.pyplot import *
import matplotlib
from scipy import linalg
# import networkx as nx
import sys
from numpy import kron
import operator
from numpy.core.umath import absolute
from basic_units import radians, degrees, cos
import math 
rcParams['figure.figsize'] = 11, 8
matplotlib.rcParams.update({'font.size': 15})

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

def plotSearch(N,theta,probT,steps,configVec):
    stepAux = []
    stepList = []
    plotName = ""
    for step in steps:
        for i in range(1,step+1):
            stepAux.append(i)
        stepList.append(stepAux)
        stepAux = []
    for steps,walk,config,n in zip(stepList,probT,configVec,N):
        plot(walk,color=config[0],linestyle=config[1],label="N=%s"%(n))
        vlines(max(steps),0,walk[-1],color=config[0],linestyle=config[2])
        legend()
        xlabel("Number of steps")
        ylabel("Probability of the marked element")
    for n in N:
        plotName+=str(n)
    savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Simulations/StagQuantumWalk/Search/'+plotName)
    clf()

def plotTheta(N,theta,probT):
    thetaDict = {}
    tList = []
    tList2 = []
    for prob,th in zip(probT,theta):
        #print("Max prob:%s \t\t Theta:%s "%(max(prob),th))
        tList.append((th,max(prob)))
        
    x,y = zip(*tList)
    xlim(0,np.pi)
    xticks(np.linspace(0, np.pi, 5),['0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'])
    xlabel("Value of θ")
    ylabel("Maximum probability of the marked element for N=%s"%N[0])
    vlines(np.pi/2,0,max(y),linestyle = '--', color = 'b')
    plot(x,y,xunits=radians)
    show()
    return tList

def plotTheta2(N,theta,probT):
    thetaDict = {}
    tList = []
    tList2 = []
    for prob,th in zip(probT,theta):
        tList.append((th,max(prob)))
    return tList


def plotMultipleThetas(N,configVec):
    probTMTheta = []
    probT = []
    thetaDist = []
    plotName=""
    for n in N:
        tVec = spaceGen(n)
        H = completeTessList(n)
        oracle = oracleList(n,0)
        theta = np.linspace(0,np.pi,Samples).tolist()
        # theta = np.linspace(0,2*np.pi,Samples).tolist()
        probT= staggeredSearchList(n, tVec, 0, oracle, H, theta,configVec)
        thetaDist = plotTheta2(n,theta,probT)
        probTMTheta.append(thetaDist)
        probT = []
        thetaDist = []
        plotName+=str(n[0])
    for probTheta,n,config in zip(probTMTheta,N,configVec):
        x,y = zip(*probTheta)
        plot(x,y,color=config[0],linestyle=config[1],label="N=%s"%(n[0]))
        legend()
        xlim(0,np.pi)
        xticks(np.linspace(0, np.pi, 5),['0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'])
       # xticks(np.linspace(0, 2*np.pi, 9),['0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$','$5\pi/4$','$3\pi/2$','$7\pi/2$','$2\pi$'])
        xlabel("Value of θ")
        ylabel("Maximum probability of the marked element")
        vlines(np.pi/2,0,max(y),linestyle = '--', color = 'b')
       # vlines([np.pi/2,3*np.pi/2],0,max(y),linestyle = '--', color = 'b')
    savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Simulations/StagQuantumWalk/Search/Theta'+plotName)
    clf()

    
def spaceGen(N):
    stepVec = []
    tVec = []
    for n in N:
        idealSteps = floor((pi/4)*sqrt(n))
        stepVec.append(int(idealSteps))
    return stepVec

def staggeredSearchList(N,tSpace,marked,oracleList,completeTessList,thetas,configVec):
    prob = []
    probT = []
    for n,oracle,tess,steps,theta in zip(N,oracleList,completeTessList,tSpace,thetas):
        evol = evo(theta,tess,oracle)
        psiN = init(n)
        prob += [absolute(psiN[marked][0])**2]
        for step in range(1,steps+1):
            psiN = evol.dot(psiN)
            prob += [(absolute(psiN[marked][0])**2)]
            # print(prob)
           # pairs.append(((absolute(psiN[marked][0])**2),step))
           # stepsAux.append(step)
        probT.append(prob)
        #stepsAuxT.append(stepsAux)
        #print("Experimental steps:%s\tTheoretical Steps:%s\n"%(max(pairs),(pi/4)*sqrt(n)))
    #    pairs = []
    #    stepsAux = []
        prob = []
    return probT

def staggeredSearch(N,U,steps,marked):
    psiN = init(N)
    probs = zeros((steps+1,1))
    for t in range(1,steps+1):
        psiN = U.dot(psiN)
        probAux = ampToProb(N,psiN,marked)
        probs[t] = probAux[marked]
    return psiN,probs

Samples = 200
plotThetaN = [256]*Samples

plotSearchN=[16,32,64]

plotMultipleThetaN = [[16]*Samples,[32]*Samples,[64]*Samples]
# print(plotMultipleThetaN)

plotSearchtheta = [(pi/2)]*len(plotSearchN)
marked = 0
#plotThetatheta =np.linspace(0,np.pi,Samples).tolist()

#tVecTheta = spaceGen(plotThetaN)
tVecSearch = spaceGen(plotSearchN)

#HTheta=completeTessList(plotThetaN)
HSearch=completeTessList(plotSearchN)

#oracleTheta = oracleList(plotThetaN,marked)
oracleSearch = oracleList(plotSearchN,marked)


colors = ['r','b','g','k']
lines = ['-','-','-' ,'-']
lines2 = ['--','--','--','--']
configVec = zip(colors,lines,lines2)

# probTTheta= staggeredSearchList(plotThetaN, tVecTheta, marked, oracleTheta, HTheta, plotThetatheta,configVec)

probTSearch = staggeredSearchList(plotSearchN, tVecSearch, marked, oracleSearch, HSearch, plotSearchtheta,configVec)
#plotSearch(plotSearchN,plotSearchtheta,probTSearch,tVecSearch,configVec)

plotMultipleThetas(plotMultipleThetaN,configVec)
#plotMTN2 = [[16]*100,[32]*100,[64]*100]

#plotMultipleThetas(plotMTN2,configVec)
