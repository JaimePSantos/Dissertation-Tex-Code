from scipy import linalg
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from IPython.display import HTML, display, Image
from math import (log,ceil,floor)
from scipy.fft import fft, ifft
from scipy.linalg import dft, inv, expm, norm
from numpy.linalg import matrix_power
from matplotlib import animation
from matplotlib import pyplot as plt

def circulantAdjacency(n,v):
    iv = list(range(0,n))
    av = list(range(0,n-1))
    C = np.zeros([n,n])
    for z in range(n):
        C[z,0] = v[iv[z]]    
    for x in range(1,n):
        av = iv[0:-1]
        iv[0] = iv[-1]
        iv[1::] = av
        for y in range(0,n):
            C[y,x] = v[iv[y]]
    return C

def initStateCont(N,initCond):
    psi0 = np.zeros((N,1))
    for x in initCond:
        psi0[x] = 1 / len(initCond)
    return psi0

def final_state(Op,psi0):
    psiN = Op.dot(psi0)
    return psiN
def ct_evo(H,t,gamma):
    U = linalg.expm(-1j*gamma*H*t)
    return U

def prob_vec(psiN,N):
    probs = np.zeros((N,1))
    for x in range(N):
        probs[x]=psiN[x]*np.conjugate(psiN[x]) 
    return probs

def ctqwalk(N, A, t, gamma, initState):
    psi0 = initState
    U = ct_evo(A,t,gamma)
    psiN = final_state(U,psi0)
    probvec = prob_vec(psiN,N)
    return probvec

def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    N = 200
    cCycle = [0,1] + [0 for x in range(N-3)] + [1]
    denom = 2
    gamma = 1/(denom*np.sqrt(2))
    A = circulantAdjacency(N,cCycle)
    initCond = [int(N/2),int(N/2)+1]
    psi0 = initStateCont(N,initCond)
    U = ct_evo(A,i,gamma)
    psiN = final_state(U,psi0)

    probs = np.zeros((N,1))
    for k in range(N):
        probs[k]=psiN[k]*np.conjugate(psiN[k])
    x = np.arange(-N/2,N/2)
    y = probs
    line.set_data(x, y)
    return line,

fig = plt.figure()
ax = plt.axes(xlim=(-100, 100), ylim=(0, 0.03))
line, = ax.plot([], [], lw=2)
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=500, interval=200, blit=True)
anim.save('ContQuantumWalk.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

