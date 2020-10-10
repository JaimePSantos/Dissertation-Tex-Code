from numpy import *
from matplotlib.pyplot import *
from graphs import *


def coins(Matrix):
    if Matrix == "H":
        coin = array([[1/sqrt(2), 1/sqrt(2)], [1/sqrt(2), -1/sqrt(2)]])
    elif Matrix == "X":
        coin = array([0, 1], [1, 0])
    return coin


def init_state(N, Pos, Amp, CoinState):
    initstate = zeros((N, 1))
    for Node in range(0, Pos.size):
        initstate[Pos[Node], 0] = Amp[Node]
    g = kron(CoinState, initstate)
    return g


def init_cond(init, state0, state1):
    if init == "0":
        psi0 = array([[1], [0]])
    if init == "1":
        psi0 = array([[0], [1]])
    if init == "01":
        psi0 = array([[1/sqrt(2)], [(1*1j)/sqrt(2)]])
    return psi0


def walk_op(N, graph,state0,state1):
    c00= outer(state0,state0)
    c01= outer(state0,state1)
    c10= outer(state1,state0)
    c11= outer(state1,state1)
    Shift = graph.adjacency_matrix()
    # Shift = kron(c00,Shift1) + kron(c11,Shift1)
    print(Shift)
    return Shift


def CU_op(coin, shift, N):
    # print(coin)
    U = shift.dot(kron(coin, eye(4)))
    return U


def final_state(U, psi0, steps):
    for t in range(0, steps):
        psi0 = U.dot(psi0)
    return psi0


def prob_vec(psiN, N):
    probs = zeros((N, 1))
    for x in range(N):
        probs[x] = psiN[x]*conjugate(psiN[x]) + \
            psiN[N+x]*conjugate(psiN[N+x])  # duvida aqui
    return probs


def plotqw(N, prob):
    matplotlib.rcParams.update({'font.size': 14})
    fig = figure()
    ax = fig.add_subplot(1, 1, 1)  # add_subplot(nrows, ncols, index, **kwargs)
    x = range(N)
    plot(x, prob)  # plot the lines
    xlabel("Graph Node")
    ylabel("Probability")
    show()


def cqwalk(N, Steps, state0, state1):
    P = int((N+1)/2)
    g = Graph({})
    graph = g.compBipartite(N)
    Coin = coins("H")

    shift = walk_op(N, graph,state0,state1)

    U = CU_op(Coin, shift, N)

    coinstate = init_cond("01", state0, state1)

    amp = array([1])
    psi0 = init_state(N, array([P]), amp, coinstate)
    psiN = final_state(U, psi0, Steps)
    probvec = prob_vec(psiN, N)

    plotqw(N, probvec)


N = 8
Steps = 1

# circuit
state0 = array([1, 0])
state1 = array([0, 1])
# circuit

cqwalk(N, Steps, state0, state1)
