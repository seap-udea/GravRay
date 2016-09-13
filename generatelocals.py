from gravray import *

###################################################
#USAGE
###################################################
usage="""Generate local conditions file for a GRT analysis.

python generateinitialconditions.py <file.dirs> <nvel> <method> <qapex>,<file.vels>

Where:

  <file.dirs>: file with random directions to be used as Azimuths
               (longitude) and latitudes.  It can be generated using
               generatedirections.py script. There are a set of
               precalculated directions that could be used here and
               that are in the util/data directory

  <nvel>: number of velocities to generate per direction

  <method>: generation method, either random (following PDF) or
            regular (spaced proportionally to probability)

  <qapex>,<file.vels>: list of files with cumulative distributions for
                       velocities at different apex directions.  qapex
                       is the initial value of the qapex interval (-1
                       if it is the whole interval).

Return:

  Files:

     locals.dat: output file cointaining the initial conditions of
                 the simulation (a matrix with Az, h and velocities).

  Figures: all figures are generated at the scratch directory.
     
     initial-directions.png: 3d map of the initial directions.


Example:

  python generatelocals.py util/data/directions-r3.00e+01.data 10 regular -1,util/data/vdistrib-qapex_0_180.data 0,util/data/vdistrib-qapex_0_60.data 60,util/data/vdistrib-qapex_60_120.data 120,util/data/vdistrib-qapex_120_180.data

Output example: locals.dat


    #1:Elevation(deg)   2:Azimuth(deg)      3:v(qapex=-1)       4:v(qapex=0)        5:v(qapex=60)       6:v(qapex=120)      
    +6.6036e+00         +1.9436e+02         +1.2864e+01         +1.1968e+01         +1.3488e+01         +1.2249e+01         
    +6.6036e+00         +1.9436e+02         +1.3983e+01         +1.2694e+01         +1.4824e+01         +1.2804e+01         
    +6.6036e+00         +1.9436e+02         +1.5054e+01         +1.3461e+01         +1.6042e+01         +1.3233e+01         

"""

###################################################
#INPUTS
###################################################
qshow=0
iarg=1
try:
    randir=argv[iarg];iarg+=1
    nvel=int(argv[iarg]);iarg+=1
    method=argv[iarg];iarg+=1

    qapexs=[]
    velcums=[]
    for velopt in argv[iarg:]:
        qapex,velcum=velopt.split(",")
        qapexs+=[qapex]
        velcums+=[velcum]
    napex=len(qapexs)
except:
    print usage
    exit(1)

###################################################
#OPEN FILE
###################################################
f=open("locals.dat","w")
f.write("%-20s%-20s"%("#1:Elevation(deg)","2:Azimuth(deg)"))
c=3
for i in xrange(napex):
    s="%d:v(qapex=%s)"%(c,qapexs[i])
    f.write("%-20s"%s)
    c+=1
s="%d:End"%c
f.write("%-20s"%s)
f.write("\n")
print "Number of apex velocities:",napex

###################################################
#READ DIRECTIONS
###################################################
datadir=np.loadtxt(randir)
bs=datadir[:,3]
ndirs=datadir[bs>0].shape[0]
print "Number of random directions:",ndirs

###################################################
#GENERATE VELOCITIES
###################################################
vels=np.zeros((nvel,napex))
for j in xrange(napex):
    source=velcums[j]

    #RANDOM 
    if method=="random":
        vs=generateVelocities(source,nvel)

    #UNIFORMLY DISTRIBUTED IN CUMMULATIVE 
    else:
        du=1./nvel
        u=np.linspace(du,0.999,nvel)
        velcum=np.loadtxt(source)
        ifvelcum=interp1d(velcum[:,1],velcum[:,0])
        vs=ifvelcum(u)

    vels[:,j]=vs
print "Number of velocities:",nvel

###################################################
#GENERATE INITIAL CONDITIONS
###################################################
n=0
Azs=[]
Els=[]
for direction in datadir[:,2:]:
    l=np.mod(direction[0]*RAD,360)
    b=direction[1]*RAD
    if b<0:continue
    Azs+=[l]
    Els+=[b]

    for i in xrange(nvel):
        f.write("%-+20.4e%-+20.4e"%(b,l))
        for j in xrange(napex):
            f.write("%-+20.4e"%(vels[i,j]))
        f.write("%-+20d\n"%(-1))
        n+=1
    #f.write("#\n")

print "%d initial conditions generated..."%n


###################################################
#PLOT RANDOM DIRECTION IN 3D
###################################################
plt.close("all")

Els=np.array(Els)
Azs=np.array(Azs)

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
