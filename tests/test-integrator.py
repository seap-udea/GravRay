#!/usr/bin/env python
from gravray import *
qcalc=1
qcomp=1

try:
    qcalc=int(argv[1])
    qcomp=int(argv[2])
except:
    pass

#INITIAL DATE
inidate="07/19/2015 00:00:00.000 UTC"
T=-1.0 #Years
nsteps=10 #Number of steps
dt=(T/nsteps)*YEAR #Output step size

#INITIAL POSITION
print "Computing initial position..."

#"""
#To use the moon inactivate the MOON in the objects.cpp file
out=System("rm whereisit.exe ; make whereisit.exe")
print out
cmd="./whereisit.exe MOON '%s' 2> /dev/null"%inidate
#"""

"""
#To run an asteroid uncomment the moon
out=System("rm whereisthisasteroid.exe ; make whereisthisasteroid.exe")
print out
cmd="./whereisthisasteroid.exe EROS '%s' 2> /dev/null"%inidate
#"""

print "Running:",cmd
inistate=System(cmd);

state=out2state(inistate)
tini=state[0]

if qcalc:
    #NUMERICAL INTEGRATION
    print "Integrating orbit..."
    out=System("rm wherewillitbe.exe ; make wherewillitbe.exe")
    print out
    cmd="./wherewillitbe.exe %s %f %d 2> /dev/null"%(inistate,T,nsteps)
    print "Running:",cmd
    out=System(cmd)

if qcomp:
    #INTERMEDIATE POSITIONS
    print "Checking integration..."
    data=np.loadtxt("ray.dat")
    dstate=[]
    for i in xrange(data.shape[0]):
        print "i = ",i
        tref=data[i,0]*YEAR
        state_integ=data[i,1:7]*STATECONV

        #INSTANTANEOUS TIME
        t=tini+tref
        print "\t","Computing position at ET = %.17e"%t

        #COMPUTE SPICE EPHEMERIS
        out=System("./whereisit.exe MOON ET %.9f 2> /dev/null"%t);
        # out=System("./whereisthisasteroid.exe EROS ET %.9f 2> /dev/null"%t);
        state_spice=out2state(out,ini=1)
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
