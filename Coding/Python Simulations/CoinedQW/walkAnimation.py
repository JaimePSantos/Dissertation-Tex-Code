import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

def coins(Matrix):
    if Matrix == "H":
        coin = np.array([[1/np.sqrt(2) , 1/np.sqrt(2)],[1/np.sqrt(2) , -1/np.sqrt(2)]])
    elif Matrix == "X":
        coin = np.array([0,1],[1,0])
    return coin

def init_state(N,Pos,Amp,CoinState):
    initstate = np.zeros((N,1))
    for Node in range(0,Pos.size):
        initstate[Pos[Node],0] = Amp[Node]
    g = np.kron(CoinState,initstate)
    return g

def init_cond(init,state0,state1):
    if init == "0":
        psi0 = np.array([[1],[0]])
    if init == "1":
        psi0 = np.array([[0],[1]])
    if init == "01":
        psi0 = np.array([[1/np.sqrt(2)],[(-1*1j)/np.sqrt(2)]])
    return psi0

def walk_op(N,state0,state1):
    c00= np.outer(state0,state0)
    c01= np.outer(state0,state1)
    c10= np.outer(state1,state0)
    c11= np.outer(state1,state1)
    ShiftPlus = np.roll(np.eye(N),1,axis=0)
    ShiftMinus = np.roll(np.eye(N),-1,axis=0)
    Shift = np.kron(c00,ShiftPlus) + np.kron(c11,ShiftMinus) 
    return Shift

def CU_op(coin, shift,N):
    U = shift.dot(np.kron(coin,np.eye(N)))
    return U

def final_state(U,psi0,steps):
    for t in range(0,steps):
        psi0=U.dot(psi0)
    return psi0

def prob_vec(psiN,N):
    probs = np.zeros((N,1))
    for x in range(N):
        probs[x]=psiN[x]*np.conjugate(psiN[x]) + psiN[N+x]*np.conjugate(psiN[N+x]) #duvida aqui
    return probs

fig = plt.figure()
ax = plt.axes(xlim=(-100, 100), ylim=(0, 0.2))
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    N = 200
    P = int((N+1)/2)
    Coin = coins("H")
    state0 = np.array([1,0])
    state1 = np.array([0,1])
    shift= walk_op(N,state0,state1)
    init = '1'
    initcond = init_cond(init,state0,state1)
    U = CU_op(Coin,shift,N)
    coinstate = initcond
    amp = np.array([1])
    psi0 = init_state(N,np.array([P]),amp,coinstate)
    psiN = final_state(U,psi0,i)
    probs = np.zeros((N,1))
    for k in range(N):
        probs[k]=psiN[k]*np.conjugate(psiN[k]) + psiN[N+k]*np.conjugate(psiN[N+k])
    x = np.arange(-N/2,N/2)
    y = probs
    line.set_data(x, y)
    return line,

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=500, interval=200, blit=True)
anim.save('QuantumWalkAsym.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
print("bla")
#plt.show()
