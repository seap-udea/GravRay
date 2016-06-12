from gravray import *

###################################################
#USAGE
###################################################
usage="""Generate initial conditions file for a GRT analysis.

python generateinitialconditions.py <file.dirs> <file.vels> <file.initial>

Where:

  <file.dirs>: file with random directions to be used as Azimuths
  (longitude) and latitudes.  It can be generated using
  generatedirections.py script. There are a set of precalculated
  directions that could be used here and that are in the util/data
  directory

  <file.vels>: file with the list of velocities that will be used in
  the simulation.

  <file.initial>: output file cointaining the initial conditions of
  the simulation (a matrix with Az, h and v).

"""

###################################################
#INPUTS
###################################################
iarg=1
try:
    randir=argv[iarg];iarg+=1
    ranvel=argv[iarg];iarg+=1
    inifile=argv[iarg];iarg+=1
except:
    print usage
    exit(1)
try:
    qshow=int(argv[iarg])
    iarg+=1
except:qshow=False

f=open(inifile,"w")
datadir=np.loadtxt(randir)
datavel=np.loadtxt(ranvel)

print "Number of random directions:",datadir.shape[0]
print "Number of random velocities:",datavel.shape[0]

###################################################
#GENERATE INITIAL CONDITIONS
###################################################
f.write("%-20s%-20s%-20s\n"%("#1:Elevation(deg)",
                             "2:Azimuth(deg)",
                             "3:Velocity(km/s)"))

n=0
Azs=[]
Els=[]
for direction in datadir[:,2:]:
    l=np.mod(direction[0]*RAD,360)
    b=direction[1]*RAD
    if b<0:continue
    Azs+=[l]
    Els+=[b]
    for v in datavel:
        f.write("%-+20.4e%-+20.4e%-+20.4e\n"%(b,l,v))
        n+=1
print "%d initial conditions generated..."%n
Els=np.array(Els)
Azs=np.array(Azs)

###################################################
#PLOT RANDOM DIRECTION IN 3D
###################################################
plt.close("all")
print "3D Map of directions..."
xs=np.cos(Els*DEG)*np.cos(Azs*DEG)
ys=np.cos(Els*DEG)*np.sin(Azs*DEG)
zs=np.sin(Els*DEG)

fig=plt.figure()
ax3d=plt3d(fig)
ax3d.set_aspect('equal')

u,v=np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x=np.cos(u)*np.sin(v)
y=np.sin(u)*np.sin(v)
z=np.cos(v)
ax3d.plot_wireframe(x, y, z, color="k")
ax3d.plot(xs,ys,zs,'o')
ax3d.plot_surface(x,y,z,rstride=1,cstride=1,color='c',alpha=1,linewidth=0)
plt.savefig("scratch/initial-directions.png")
if qshow:plt.show()

###################################################
#END GENERATION
###################################################
f.close()
