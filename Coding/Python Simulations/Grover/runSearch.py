#GroverSearch

import numpy as np
from matplotlib.pyplot import *
import matplotlib
from scipy import linalg
rcParams['figure.figsize'] = 11, 8
matplotlib.rcParams.update({'font.size': 15})

def init(N):
    psi0 = np.ones((N,1))/ np.sqrt(N)
    return psi0

def oracle(N,marked):
    oracle = np.eye(N)
    for mark in marked:
        oracle[mark][mark] = -1
    return oracle

def diffusion(N):
    ketSuper = np.ones((N,1))/ np.sqrt(N)
    diff = 2 * np.outer(ketSuper,ketSuper) - np.eye(N)
    return diff

def unitary(N,marked):
    orac = oracle(N,marked)
    diff = diffusion(N)
    return np.dot(diff,orac)

def spaceGen(N):
    stepVec = []
    tVec = []
    for n in N:
        idealSteps =np.floor((np.pi/4)*np.sqrt(n))
        stepVec.append(int(idealSteps))
    return stepVec

def spaceGenMultiple(N,numberOfMarked):
    stepVec = []
    tVec = []
    for n in N:
        idealSteps =np.floor((np.pi/4)*np.sqrt(n/numberOfMarked))
        stepVec.append(int(idealSteps))
    return stepVec

def groverSearch(N,stepSpace,marked):
    prob = []
    probT = []
    for n,steps in zip(N,stepSpace):
        u = unitary(n,marked)
#        print("N=%s\tSteps=%s"%(n,steps))
        psiN=init(n)
        prob += [np.absolute(psiN[marked][0][0])**2]
        for step in range(1,steps+1):
            psiN = np.dot(u,psiN)
            prob+=[np.absolute(psiN[marked][0][0]**2)]
        probT.append(prob)
        prob = []
    return probT

def plotSearch(N,probT,steps,configVec):
    stepAux = []
    plotName = ""
    stepList = []
    for step in steps:
        for i in range(1,step+1):
            stepAux.append(i)
        stepList.append(stepAux)
        stepAux = []
    for n,steps,config,walk in zip(N,stepList,configVec,probT):
        plot(walk,color=config[0],linestyle=config[1],label="N=%s"%(n))
        vlines(max(steps),0,walk[-1],color=config[0],linestyle=config[2])
        legend()
        xlabel("Number of steps")
        ylabel("Probability of the marked element")
    for n in N:
        plotName+=str(n)
    savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Simulations/Grover/GroverOneMarked'+plotName)
    clf()

def groverSearchMultiple(N,stepSpace,marked):
    prob = 0
    probT = []
    probTAux = []
    i =0
    for n,steps in zip(N,stepSpace):
        u = unitary(n,marked)
#        print("MultipleMarked N=%s\tSteps=%s"%(n,steps))
        psiN=init(n)
        prob=(np.absolute((psiN[marked][0][0])**2)*len(marked))
        probTAux.append(prob)
        for step in range(1,steps+1):
            psiN = np.dot(u,psiN)
            for mark in marked: 
                prob+=np.absolute((psiN[mark][0]))**2
#            print("inside "+str(prob))
#            print()
            probTAux.append(prob)
            prob=0
#        print("outside"+str(psiN))
        probT.append(probTAux)
        probTAux = []
#    print(probT)
    return probT

def plotSearchMultiple(N,probT,steps,configVec):
    stepAux = []
    plotName = ""
    stepList = []
    for step in steps:
        for i in range(1,step+1):
            stepAux.append(i)
        stepList.append(stepAux)
        stepAux = []
    for n,steps,config,walk in zip(N,stepList,configVec,probT):
        plot(walk,color=config[0],linestyle=config[1],label="N=%s"%(n))
        vlines(max(steps),0,walk[-1],color=config[0],linestyle=config[2])
        legend()
        xlabel("Number of steps")
        ylabel("Total probability of the marked elements")
    for n in N:
        plotName+=str(n)
    savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Simulations/Grover/GroverMultipleMarked'+plotName)
    clf()

def markedList(N):
    elementListAux1 = []
    elementListAux2 = []
    elementList = []
    for n in N:
        markedElementRange = np.arange(int(n/4)+1).tolist()
        for element in markedElementRange:
            for i in range(element):
                elementListAux1.append(i)
            if len(elementListAux1)>0:
                elementListAux2.append(elementListAux1)
            elementListAux1 =[]
        # elementListAux2.pop(0)
        elementList.append(elementListAux2)
        elementListAux2 = []
    return elementList

def singleShotGrover(N,markedListListList):
    prob = 0
    probTAux = []
    probT = []
    for n,markedListList in zip(N,markedListListList):
        psiN = np.zeros(n)
        for markedList in markedListList:
            u = unitary(n,markedList)
            psi0=init(n)
            psiN = np.dot(u,psi0)
            for marked in markedList:
                prob+=np.absolute(psiN[marked][0])**2
            probTAux.append(prob)
            prob=0
        probT.append(probTAux)
        probTAux = []
    return probT

def plotSingShot(N,walkList,configVec):
    plotName=""
    for walk,config,n in zip(walkList,configVec,N):
        plot(np.append(np.roll(walk,1),walk[int(n/4)-1]),color=config[0],linestyle=config[1],label="N=%s"%(n))
        # plot(walk,color=config[0],linestyle=config[1],label="N=%s"%(n))
        xlim(1,int(n/4)+1)
        dim = np.arange(1,int(n/4)+1,3)
        xticks(dim)
        vlines(int(n/4),0,walk[-1],color=config[0],linestyle=config[2])
        legend()
        xlabel("Number of marked elements")
        ylabel("Total probability of marked elements")
    for n in N:
       plotName+=str(n) 
    savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Simulations/Grover/GroverSingleShot'+plotName)
    clf()


colors = ['r','b','g','k']
lines = ['-','-','-' ,'-']
lines2 = ['--','--','--','--']
configVec = zip(colors,lines,lines2)

NSingShot =[64,128,256]
markedSingShot = markedList(NSingShot)
groverSingle = singleShotGrover(NSingShot,markedSingShot)
plotSingShot(NSingShot,groverSingle,configVec)

configVec2 = zip(colors,lines,lines2)
N=[32,64,128,256]
marked=[0]
steps=spaceGen(N)
grover = groverSearch(N,steps,marked)
plotSearch(N,grover,steps,configVec2)

NMMarked = [32,64,128,256]
mMarked = [0,1]
configVec3 = zip(colors,lines,lines2)
stepsMMarked = spaceGenMultiple(NMMarked,len(mMarked))
groverMMarked = groverSearchMultiple(NMMarked,stepsMMarked,mMarked)
plotSearchMultiple(NMMarked,groverMMarked,stepsMMarked,configVec3)

