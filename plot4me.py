#!/usr/bin/env python
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D as plt3d
import numpy as np

"""
Fields: 
1:t
2:x
3:y
4:z
5:vx
6:vy
7:vz
8:q
9:e
10:i
11:W
12:w
13:M
14:t0
15:mu
"""                     

def plotRay():
    data=np.loadtxt("ray.dat")

    fig=plt.figure()
    ax=fig.gca()
    ax.plot(data[:,0],data[:,7])
    fig.savefig("plots/elements.png")
    plt.close("all")

    fig3d=plt.figure()
    ax=plt3d(fig3d)
    ax.plot(data[:,1],data[:,2],data[:,3])
    fig3d.savefig("plots/trajectory3d.png")
    plt.show()

plotRay()
