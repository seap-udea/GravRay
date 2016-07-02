#!/usr/bin/env python
from gravray import *

#DO YOU WANT TO PROCEED
qcalc=1

#WHICH OBJECT
objtype=argv[1]
objname=argv[2]

if objtype=="Planet":
    posprogram="whereisit"
else:
    posprogram="whereisthisasteroid"

try:
    qcalc=int(argv[3])
except:
    pass

#############################################################
#INITIAL DATE
#############################################################
inidate="07/19/2015 00:00:00.000 UTC"

#############################################################
#INTEGRATION PROPERTIES
#############################################################
T=+5.0 #Years
nsteps=100 #Number of steps
dt=(T/nsteps)*YEAR #Output step size

if qcalc:

    #############################################################
    #INITIAL POSITION
    #############################################################
    print "Computing initial position..."

    out=System("make %s.exe"%(posprogram))
    print "Output:\n",out
    cmd="./%s.exe %s '%s' > /dev/null"%(posprogram,objname,inidate)
    print "Running:",cmd
    out=System(cmd);
    dst=out2dict(out)
    state=dst["STATE(6)"]
    tini=dst["TDB"]
    
    #############################################################
    #NUMERICAL INTEGRATION
    #############################################################
    print "Integrating orbit forward..."
    out=System("make wherewillitbe.exe")
    print "Output:\n",out
    cmd="./wherewillitbe.exe %.9e %s %.17e %d &> /tmp/err.log"%(tini,vec2str(state),T,nsteps)
    print "Running:",cmd
    out=System(cmd)
    System("cp ray.dat scratch/ray-forward.dat")
    
    #############################################################
    #PROPAGATING BACKWARDS
    #############################################################
    data=np.loadtxt("ray.dat")
    state=data[-1,1:7]
    tend=data[-1,0]

    cmd="./wherewillitbe.exe %.9e %s %.17e %d &> /tmp/err.log"%(tend,vec2str(state),-T,nsteps)
    print "Running:",cmd
    out=System(cmd)
    System("cp ray.dat scratch/ray-backward.dat")

data_forward=np.loadtxt("scratch/ray-forward.dat")
data_backward=np.loadtxt("scratch/ray-backward.dat")
data_dif=np.abs(data_forward[:,1:]-data_backward[::-1,1:])/np.abs(data_forward[:,1:]+data_backward[::-1,1:])/2
tini=data_forward[0,0]
ts=(data_forward[:,0]-tini)/DAY

fig=plt.figure()
ax=fig.gca()
ax.plot(ts,data_dif[:,0],label='x')
ax.plot(ts,data_dif[:,1],label='y')
ax.plot(ts,data_dif[:,2],label='z')
ax.legend(loc='best')
ax.set_yscale("log")
ax.set_title("Position, T = %.2f years"%T)
fig.savefig("scratch/dray-pos-direction.png")

fig=plt.figure()
ax=fig.gca()
ax.plot(ts,data_dif[:,3],label='x')
ax.plot(ts,data_dif[:,4],label='y')
ax.plot(ts,data_dif[:,5],label='z')
ax.legend(loc='best')
ax.set_yscale("log")
ax.set_title("Velocity, T = %.2f years"%T)
fig.savefig("scratch/dray-vel-direction.png")
