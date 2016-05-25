#!/usr/bin/env python
from gravray import *

#DO YOU WANT TO INTEGRATE?
qcalc=1

#DO YOU WANT TO COMPARE INTEGRATION WITH SPICE
qcomp=1

#WHICH OBJECT
objtype=argv[1]
objname=argv[2]

if objtype=="Planet":
    posprogram="whereisit"
else:
    posprogram="whereisthisasteroid"

#############################################################
#GET INPUT
#############################################################
try:
    qcalc=int(argv[3])
    qcomp=int(argv[4])
except:
    pass

#############################################################
#INITIAL DATE
#############################################################
inidate="07/19/2015 00:00:00.000 UTC"

#############################################################
#INTEGRATION PARAMETERS
#############################################################
T=+1.0 #Years
nsteps=10 #Number of steps
dt=(T/nsteps)*YEAR #Output step size


if qcalc:
    #############################################################
    #INITIAL POSITION
    #############################################################
    print "Computing initial position..."

    #To use the moon inactivate the MOON in the objects.cpp file
    out=System("make %s.exe"%(posprogram))
    print "Output:\n",out
    cmd="./%s.exe %s '%s' > /dev/null"%(posprogram,objname,inidate)
    print "Running:",cmd
    out=System(cmd);
    #print "Output:\n",out
    dst=out2dict(out)

    state=dst["STATE(6)"]
    tini=dst["TDB"]

    #############################################################
    #NUMERICAL INTEGRATION
    #############################################################
    print "Integrating orbit..."
    out=System("make wherewillitbe.exe")
    print "Output:\n",out
    cmd="./wherewillitbe.exe %.9e %s %.17e %d 2> /dev/null"%(tini,vec2str(state),T,nsteps)
    print "Running:",cmd
    out=System(cmd)
    print "Output:\n",out

if qcomp:
    #############################################################
    #SPICE POSITIONS
    #############################################################
    print "Checking integration..."
    data=np.loadtxt("ray.dat")
    dstate=[]
    for i in xrange(data.shape[0]):

        print "i = ",i
        t=data[i,0]
        state_integ=data[i,1:7]

        #========================================
        #INSTANTANEOUS TIME
        #========================================
        print "\t","Computing position at ET = %.17e"%t

        #COMPUTE SPICE EPHEMERIS
        out=System("./%s.exe %s ET %.9f > /dev/null"%(posprogram,objname,t));
        dst=out2dict(out)
        state_spice=dst["STATE(6)"]
        print "\t","SPICE:",vec2str(state_spice)

        #INTEGRATED VALUES
        print "\t","INTEG:",vec2str(state_integ)

        #CALCULATE DIFFERENCES
        d=norm(state_spice[0:3])
        v=norm(state_spice[3:])
        normstate=np.array([d,d,d,v,v,v])

        adif=np.abs((state_spice-state_integ)/normstate)
        ndif=(state_spice-state_integ)/normstate
        print "\t","DIFER:",vec2str(adif)

        dstate+=[[t]+adif.tolist()+ndif.tolist()]

    np.savetxt("scratch/dray.dat",dstate)

#PLOTTING DIFFERENCES
print "Plotting differences..."
dstate=np.loadtxt("scratch/dray.dat")
data=np.loadtxt("ray.dat")
tini=dstate[0,0]
ts=(dstate[:,0]-tini)/DAY

fig=plt.figure(figsize=(8,6))
ax=fig.gca()

ax.plot(ts,dstate[:,1],label='x')
ax.plot(ts,dstate[:,2],label='y')
ax.plot(ts,dstate[:,3],label='z')

ax.set_yscale("log")
ax.set_ylim((1e-12,1e-4))
ax.legend(loc='best')
fig.savefig("scratch/dray-pos.png")

fig=plt.figure(figsize=(8,6))
ax=fig.gca()

ax.plot(ts,dstate[:,4],label='vx')
ax.plot(ts,dstate[:,5],label='vy')
ax.plot(ts,dstate[:,6],label='vz')

ax.set_yscale("log")
ax.set_ylim((1e-12,1e-4))
ax.legend(loc='best')
fig.savefig("scratch/dray-vel.png")

fig3d=plt.figure()
ax=plt3d(fig3d)
ax.plot(dstate[:,7],dstate[:,8],dstate[:,9],'ko')
fig3d.savefig("scratch/dray-dist-pos.png")

fig3d=plt.figure()
ax=plt3d(fig3d)
ax.plot(dstate[:,10],dstate[:,11],dstate[:,12],'ko')
fig3d.savefig("scratch/dray-dist-vel.png")

#PLOTTING TRAJECTORY
plt.close("all")
fig3d=plt.figure()
ax=plt3d(fig3d)
ax.plot(data[:,1],data[:,2],data[:,3])
ext=max(np.abs(data[:,1:4].min()),data[:,1:4].max())
ax.set_xlim(-ext,ext)
ax.set_ylim(-ext,ext)
ax.set_zlim(-ext,ext)
fig3d.savefig("scratch/dray3d.png")
#plt.show()
