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
    datai=np.loadtxt("rayi.dat")

    fig=plt.figure()
    ax=fig.gca()
    ax.plot(data[:,0],data[:,7])
    fig.savefig("scratch/elements.png")
    plt.close("all")

    fig3d=plt.figure()
    ax=plt3d(fig3d)
    ax.plot(data[:,1],data[:,2],data[:,3])
    ax.plot(datai[:,1],datai[:,2],datai[:,3],'r-')
    fig3d.savefig("scratch/trajectory3d.png")

    ext=max(np.abs(data[:,1:4].min()),data[:,1:4].max())
    ax.set_xlim(-ext,ext)
    ax.set_ylim(-ext,ext)
    ax.set_zlim(-ext,ext)
    plt.show()

plotRay()
