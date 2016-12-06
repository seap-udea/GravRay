from gravray import *

###################################################
#USAGE
###################################################
usage="""Generate random impact velocities for a GRT analysis.

   python generatevelocities.py <nvel> <method> [<qapex1>,<file1.vels> <qapex2>,<file2.vels> ...]

Where:

  <nvel>: number of velocities to generate per direction

  <method>: generation method, either random (following PDF),
            regular (spaced proportionally to probability) or uniform.

  <qapex>,<file.vels>: list of files with cumulative distributions for
                       velocities at different apex directions.  qapex
                       is the initial value of the qapex interval (-1
                       if it is the whole interval).

Return:

  Files:

     velocities.dat: matrix with velocities

Example:

  python generatevelocities.py 10 regular -1,util/data/vdistrib-qapex_0_180.data

  python generatevelocities.py 10 regular -1,util/data/vdistrib-qapex_0_180.data 0,util/data/vdistrib-qapex_0_60.data 60,util/data/vdistrib-qapex_60_120.data 120,util/data/vdistrib-qapex_120_180.data

"""
###################################################
#INPUTS
###################################################
qshow=0
iarg=1
try:
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
#GENERATE VELOCITIES
###################################################
vels=np.zeros((nvel,napex))
for j in xrange(napex):
    source=velcums[j]

    #RANDOM 
    if method=="random":
        vs=generateVelocities(source,nvel)

    #UNIFORMLY DISTRIBUTED IN CUMMULATIVE 
    elif method=="regular":
        du=1./nvel
        u=np.linspace(du,0.999,nvel)
        velcum=np.loadtxt(source)
        ifvelcum=interp1d(velcum[:,1],velcum[:,0])
        vs=ifvelcum(u)

    #UNIFORMLY DISTRIBUTED
    elif method=="uniform":
        velcum=np.loadtxt(source)
        vmin=velcum[:,0].min()
        vmax=velcum[:,0].max()
        vs=np.linspace(vmin,vmax,nvel)
    else:
        print "Method %s not recognized..."
        exit(1)

    vels[:,j]=vs

print "Number of velocities:",nvel

###################################################
#SAVE VELOCITIES
###################################################
#CREATE HEADER
f=open("velocities.dat","w")
c=1
f.write("#")
for i in xrange(napex):
    s="%d:v(qapex=%s)"%(c,qapexs[i])
    f.write("%-20s"%s)
    c+=1
s="%d:End"%c
f.write("%-20s"%s)
f.write("\n")

#SAVE VELOCITIES
for i in xrange(nvel):
    for j in xrange(napex):
        f.write("%-+20.4e"%(vels[i,j]))
    f.write("%-+20d\n"%(-1))
f.close()
