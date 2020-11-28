import numpy as np
from matplotlib.pyplot import *
from scipy import stats
from scipy import stats
from scipy.stats import gaussian_kde 
import seaborn as sns
import pandas as pd
rcParams['figure.figsize'] = 11, 8
matplotlib.rcParams.update({'font.size': 15})

#N=200
#dataPoints= 72 
#aux = [0 for x in range(N+1)]
#probVec = [0 for x in range (N+1)]
#x = np.random.binomial(N,0.5,dataPoints) 
##print(x)
#for i in range(N):
#    for j in range(len(x)):
#        if i == x[j]:
#            aux[i]+=1
#
#for i in range(len(aux)):
#    probVec[i] = aux[i] / max(aux) 
#print(aux)
#auxvar = 0
##for element in aux:
##    auxvar+=element
#
##print(auxvar)
##print(probVec)
##print(probVec)
##print(x)
#plot(probVec)
#show()
# randomly walk nsteps and return the x value
# starting at x=0
# each step has zero mean and a variance of 1

# TODO: Implementar equacao 3.4 do renato portugal
def walkn(nsteps):  # random walk using a normal distribution for step sizes
    r = stats.norm.rvs(size=nsteps)  # normal distribution mean=0 variance=1
    # r is a vector values randomly generated with a normal distribution
    return sum(r)  # the sum of the entire vector

def npart_walkn(npart,nsteps):
    xvec = np.zeros(0)
    for i in range(npart):
        x = walkn(nsteps)  # a single random walk value
        xvec = np.append(xvec,x)  # append each random walk to the vector
    return xvec

def multipleWalk(npart,nsteps):
    xvecList = []
    for n in nsteps:
        xVec = npart_walkn(npart,n)
        xvecList.append(xVec)
    return xvecList

def plotMultiple(nSteps,mWalks,configVec):
    plotName=""
    for n, walk,config in zip(nSteps,mWalks,configVec):
        density = gaussian_kde(walk)
        xs = np.linspace(-100,100,500)
        density.covariance_factor = lambda : .5 
        density._compute_covariance()
        plot(xs,density(xs),color=config[0],linestyle=config[1],label="Steps=%s"%n)
        legend()
        xlabel("Graph Node")
        ylabel("Probability")
        plotName+=str(n)
    savefig(r'/home/jaime/Programming/Jaime-Santos-Dissertation/Results/Simulations/ClassicalWalk/MultClassicalWalk'+plotName)
    

nsteps = 72 # number of steps
npart = 10000# number of particles (sailors) to let walk around
xvec = npart_walkn(npart,nsteps)
density = gaussian_kde(xvec)
xs = np.linspace(-100,100,200)
density.covariance_factor = lambda : .5 
density._compute_covariance()
#plot(xs,density(xs))
#show()


colors = ['r','b','g','k']
lines = ['-','-','-','-']
lines2 = ['--','--','--','--']
configVec = zip(colors,lines,lines2)

nStepList= [72,180,450]
xVecList = multipleWalk(npart,nStepList)
plotMultiple(nStepList, xVecList,configVec)

