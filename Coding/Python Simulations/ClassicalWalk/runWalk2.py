from random import randint
import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

def classicalWalk(nExp,steps):
    positions = []
    for x in range(nExp):
        k = 0
        for y in range(steps):
            a = randint(0,1)
            k+=(-1)**a
        positions.append(k)
    return positions

def init():
    line.set_data([], [])
    return line,

def animate(i):
    ne = 20000 
    steps = i+1 
    walkList = classicalWalk(ne,steps)
    probDict = {i:walkList.count(i)/ne for i in set(walkList)}
    lists = sorted(probDict.items()) # sorted by key, return a list of tuples
    x, y = zip(*lists)
    line.set_data(x, y)
    return line,

fig = plt.figure()
ax = plt.axes(xlim=(-100, 100), ylim=(0, 0.3))
line, = ax.plot([], [], lw=2)
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=300, interval=200, blit=True)
anim.save('ClassicalWalk.mp4', fps=30, extra_args=['-vcodec', 'libx264'])
#ne = 300000
#steps = 72
#walkList = classicalWalk(ne,steps)
#probDict = {i:walkList.count(i)/ne for i in set(walkList)}
#lists = sorted(probDict.items()) # sorted by key, return a list of tuples
#x, y = zip(*lists) # unpack a list of pairs into two tuples
#plt.plot(x, y)
#plt.show()
