#!/usr/bin/env python
from gravray import *
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
    ext=0
    try:
        ext=float(argv[1])
    except:
        pass

    data=np.loadtxt("ray.dat")
    objects=getScenario("scenario.dat")
    tini=data[0,0]
    ts=(data[:,0]-tini)/YEAR

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #PLOTTING ORBITAL ELEMENTS
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig=plt.figure(figsize=(8,6))
    axq=fig.add_subplot(311,ymargin=0)
    axe=fig.add_subplot(312,ymargin=0)
    axi=fig.add_subplot(313,ymargin=0)

    axq.plot(ts,data[:,9]/UL)
    axe.plot(ts,data[:,10])
    axi.plot(ts,data[:,11])
    
    axq.set_ylabel("q (AU)")
    axe.set_ylabel("e")
    axi.set_ylabel("i (deg)")
    axi.set_xlabel(r"$t-t_{\rm impact}$ (years)")
    
    #DECORATION
    axq.set_xticklabels([])
    axe.set_xticklabels([])

    for ax in axq,axe,axi:
        ax.grid()

    fig.savefig("scratch/elements.png")
    plt.close("all")

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #PLOTTING TRAJECTORY IN 3-D
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig3d=plt.figure()
    ax=plt3d(fig3d)

    #==============================
    #RAY
    #==============================
    data[:,1:7]=np.array([data[i,1:7]/STATECONV for i in xrange(data.shape[0])])
    ax.plot(data[:,1],data[:,2],data[:,3],'k-',label='Ray')

    #==============================
    #SCENARIO
    #==============================
    for objid in objects:
        if 'ts' in objid:continue
        objdata=objects[objid]
        objdata=np.array([objdata[i,:]/STATECONV for i in xrange(objdata.shape[0])])
        cond=np.array([norm(objdata[i,0:3])<=ext for i in xrange(objdata.shape[0])])
        if len(objdata[cond,0])>0:
            ax.plot(objdata[cond,0],objdata[cond,1],objdata[cond,2],
                    marker='o',mec='none',ms=1,lw=0,label='%s'%OBJECTS[objid])

    ax.legend(loc="best")
    fig3d.savefig("scratch/trajectory3d.png")

    if ext==0:
        ext=max(np.abs(data[:,1:4].min()),data[:,1:4].max())
    ax.set_xlim(-ext,ext)
    ax.set_ylim(-ext,ext)
    ax.set_zlim(-ext,ext)

    ax.set_xlabel('x (AU)')
    ax.set_ylabel('y (AU)')
    ax.set_zlabel('z (AU)')
    plt.show()

plotRay()
