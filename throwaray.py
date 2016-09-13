from gravray import *

#############################################################
#USAGE
#############################################################
usage="""Launch an object into orbit.

python throwaray.py <lat> <lon> <alt> <elev> <azim> <vimp> <date> <tspan> [npoints]

Where:

   <lat>,<lon>: latitude and longitude in degrees.
   <alt>: altitude in meters.
   <elev>,<azim>: elevation and azimuth in degrees.
   <vimp>: impact velocity in km/s
   <date>: date in format MM/DD/CCYY HH:MM:SS.dcm UTC-L
   <tspan>: integration time (in years)
   <npoints>: sampling points

Return:

   Standard output:

       Classical elements of the asymptotic orbit

   Data files:

       ray.dat: trajectory.

   Figures: all figures are generated in the scratch directory.

       ray-elements.png: orbital elements of the ray.

       ray-orbit.png: 3d realization of the ray.
"""

#############################################################
#INPUT
#############################################################
try:
    lat=float(argv[1])
    lon=float(argv[2])
    alt=float(argv[3])
    elev=float(argv[4])
    azim=float(argv[5])
    vimp=float(argv[6])
    date=argv[7]
    tspan=float(argv[8])
except:
    print usage
    exit(1)

try:
    npoints=int(argv[9])
except:npoints=100

t=0
#############################################################
#DETERMINE INITIAL PARTICLE POSITION AND VELOCITY
#############################################################
print "Determining initial position..."
out=System("./whereami.exe %.17e %.17e %.17e %.17e %.17e %.17e '%s' > /tmp/out.log"%\
           (lat,lon,alt,
            elev,azim,vimp,
            date))
props=out2dict(out)
print TAB,"TDB = %.9e"%props["TDB"]
timeIt()

#############################################################
#INTEGRATE ORBIT
#############################################################
print "Integrating orbit..."

cmd="./wherewillitbe.exe %s %s %.17e %d > /tmp/out.log"%\
    (props["TDB"],vec2str(props["ECJ2000(6)"]),tspan,npoints)
out=System(cmd)
timeIt()

print "Computing scenario position..."
out=System("./scenarioof.exe ray.dat")
timeIt()

#############################################################
#GET FINAL ORBITAL ELEMENTS
#############################################################
data=np.loadtxt("ray.dat")
elements=data[-1,9:]
print "Orbital elements:"
print TAB,"q = ",elements[0]/UL
print TAB,"e = ",elements[1]
print TAB,"i = ",elements[2]
print TAB,"W = ",elements[3]
print TAB,"w = ",elements[4]
print TAB,"M = ",elements[5]
timeIt()

#############################################################
#PLOT ELEMENTS & TRAJECTORY
#############################################################
objects=getScenario("scenario.dat")
tini=data[0,0]
ts=(data[:,0]-tini)/YEAR

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#PLOTTING ORBITAL ELEMENTS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print "Plotting elements..."
fig=plt.figure(figsize=(8,6))
axq=fig.add_subplot(311,ymargin=0)
axe=fig.add_subplot(312,ymargin=0)
axi=fig.add_subplot(313,ymargin=0)

axq.plot(ts,data[:,9]/UL)
axe.plot(ts,data[:,10])
axi.plot(ts,data[:,11])

axq.set_ylabel("$q$ (AU)")
axe.set_ylabel("$e$")
axi.set_ylabel("$i$ (deg)")
axi.set_xlabel(r"$t-t_{\rm impact}$ (years)")

#DECORATION
axq.set_xticklabels([])
axe.set_xticklabels([])

for ax in axq,axe,axi:
    ax.grid()

fig.savefig("scratch/ray-elements.png")
plt.close("all")
timeIt()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#PLOTTING TRAJECTORY IN 3-D
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print "Plotting trajectory..."
fig3d=plt.figure()
ax=plt3d(fig3d)
data[:,1:7]=np.array([data[i,1:7]/STATECONV for i in xrange(data.shape[0])])

#==============================
#SCENARIO
#==============================
ext=0
if ext==0:
    ext=0.8*max(np.abs(data[:,1:4].min()),data[:,1:4].max())

for objid in objects:
    if 'ts' in objid:continue
    objdata=objects[objid]
    objdata=np.array([objdata[i,:]/STATECONV for i in xrange(objdata.shape[0])])
    cond=np.array([norm(objdata[i,0:3])<=ext for i in xrange(objdata.shape[0])])
    if len(objdata[cond,0])>0:
        ax.plot(objdata[cond,0],objdata[cond,1],objdata[cond,2],
                '-',mec='none',ms=3,lw=1,label='%s'%OBJECTS[objid])

#==============================
#RAY
#==============================
ax.plot(data[:,1],data[:,2],data[:,3],'k--',label='Test particle')

ax.legend(loc="upper left",fontsize=10)

ax.set_xlim(-ext,ext)
ax.set_ylim(-ext,ext)
ax.set_zlim(-ext,ext)

ax.set_xlabel('x (AU)')
ax.set_ylabel('y (AU)')
ax.set_zlabel('z (AU)')
fig3d.savefig("scratch/ray-orbit.png")

timeIt()
