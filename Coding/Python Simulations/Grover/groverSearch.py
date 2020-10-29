import numpy as np
import matplotlib.pyplot as plt 
from scipy import linalg

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

def groverSearch(N,marked,steps):
    prob = []
    probT = []
    for n in N:
        u = unitary(n,marked)
        psiN=init(n)
        for t in range(steps):
            psiN = np.dot(u,psiN)
            for mark in marked:
                prob+=[np.absolute(psiN[marked][0]**2)]
            # print(prob)

        
    return prob


        
        
N=[16]
steps = 1
marked = [0,1]
N1 = [16]
steps1 = 2
N2 = [16]
steps2 = 3
N3 = [16]
steps3 = 4
# print("Oracle:\n%s\n"%oracle(N,marked))
# print("Diffusion:\n%s\n"%diffusion(N))
print("Grover evolution for %s steps and %s elements:\n%s\n"%(steps,N,groverSearch(N,marked,steps)))
print("Grover evolution for %s steps and %s elements:\n%s\n"%(steps1,N1,groverSearch(N1,marked,steps1)))
print("Grover evolution for %s steps and %s elements:\n%s\n"%(steps2,N2,groverSearch(N2,marked,steps2)))
print("Grover evolution for %s steps and %s elements:\n%s\n"%(steps2,N2,groverSearch(N3,marked,steps3)))
