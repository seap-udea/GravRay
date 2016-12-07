#!/usr/bin/env python
from gravray import *

#############################################################
#USAGE
#############################################################
usage="""Make a Parallel GRT analysis of a whole geographic region.

   python parallelgravray.py <nprocs> <date> <deg|rad> <locations_file> <deg|rad> <directions_file> <qvel> <velocities_file> <sname> [<altitude>]

Where:

   <nprocs>: number of processors.

   <date>: date of analysis.
           Use format MM/DD/CCYY HH:MM:SS.dcm UTC-L.

   <deg|rad>: indicate if directions are in degrees or radians.

   <locations_file>: file with the geographical locatins (generated
                     using eg. generatedirections.py)

   <directions_file>: file with the impacting directons (generated
                     using eg. generatedirections.py)

   <qvel>: use apex dependent velocities?

   <velocities_file>: file with initial velocities (generated
                      using eg. generatevelocities.py)

   <sname>: name of the simulation

   <altitude>: altitude where the rays start in meters (default:
               80,000 m)


"""

#############################################################
#PREPARE RUN
#############################################################
iarg=1
try:
    nprocs=int(argv[iarg]);iarg+=1

    date=argv[iarg];iarg+=1

    degloc=argv[iarg];iarg+=1
    locfile=argv[iarg];iarg+=1

    degdir=argv[iarg];iarg+=1
    dirfile=argv[iarg];iarg+=1

    qvel=int(argv[iarg]);iarg+=1
    velfile=argv[iarg];iarg+=1

    name=argv[iarg];iarg+=1
except:
    print usage
    exit(1)

#ALTITUDE
try:
    h=float(argv[iarg]);iarg+=1
except:h=8e4

#############################################################
#PREPARE
#############################################################
System("make")

#############################################################
#INPUTS ANALYSIS
#############################################################
#DATE ANALYSIS
parts=date.split(" ")
dparts=parts[0].split("/")
tparts=parts[1].split(":")
tstring="%04d%02d%02d%02d%02d%02d"%(int(dparts[2]),int(dparts[0]),int(dparts[1]),
                                    int(tparts[0]),int(tparts[1]),int(tparts[2]))
out=System("./whattimeisit.exe '%s' ET > /dev/null"%date)
t=float(out.split("\n")[4])

#MAKE STRING
dirmd5=System("md5sum %s |awk '{print $1}'"%dirfile)
velmd5=System("md5sum %s |awk '{print $1}'"%velfile)
makestr="qvel=%d & name=%s & dirmd5=%s & velmd5=%s"%(qvel,name,dirmd5,velmd5)

md5str=MD5STR(makestr,len=6)

#RUN STRING
runstr="%s-%s"%(tstring,md5str)

#OUTPUT DIRECTORY
outdir="data/grt-%s-%s"%(tstring,md5str)
System("mkdir -p %s"%outdir)
f=open(outdir+"/.config","w")
f.write("%s\n%s\n%s\n"%(md5str,tstring,makestr))
System("cp %s %s/locations.dat"%(locfile,outdir))
System("cp %s %s/directions.dat"%(dirfile,outdir))
System("cp %s %s/velocities.dat"%(velfile,outdir))
f.close()

#############################################################
#DETERMINE SIZE OF RUN
#############################################################
data=np.loadtxt(locfile)
nlocs=data.shape[0]
data=np.loadtxt(dirfile)
ndirs=data.shape[0]
data=np.loadtxt(velfile)
nvels=data.shape[0]

print "Preparing parallel run %s in %d processors..."%(runstr,nprocs)
print "\tNumber of locations: ",nlocs
print "\tNumber of directions: ",ndirs
print "\tNumber of vels: ",nvels
print "Running time estimation:"
print "\tNumber of initial conditions: ",ndirs*nvels
tploc=(ndirs*nvels)/240.0*35./60
print "\tTime per location (min): ",tploc
print "\tTotal time (hours): ",tploc*nlocs/60.
print "\tParallel time (hours): ",tploc*nlocs/60./nprocs
raw_input("Press enter to continue?..")

#############################################################
#CREATE RUN
#############################################################
locdata=np.loadtxt(locfile)
nlocs=locdata.shape[0]
nperplow=int(np.floor(nlocs/(1.0*nprocs)))
nperpup=int(np.ceil(nlocs/(1.0*nprocs)))
nperp=nperplow

iini=0
for i in xrange(nprocs):
    nproc=i+1
    iend=iini+nperp
    if i==nprocs-1:iend=nlocs
    ndata=len(locdata[iini:iend,:])
    print "\tPreparing %d locations for processor %d (i = %d - %d)..."%(ndata,nproc,iini,iend)
    np.savetxt("%s/locations-%d.dat"%(outdir,nproc),locdata[iini:iend,:])
    iini+=nperp
    if (i%2)==0:nperp=nperpup
    else:nperp=nperplow

#############################################################
#CREATE PBS SCRIPT
#############################################################
options="\"%s\" %s %s/locations-${PBS_ARRAYID}.dat %s %s/directions.dat %d %s/velocities.dat \"%s\" %e ${PBS_ARRAYID}"%(\
    date,
    degloc,outdir,
    degdir,outdir,
    qvel,outdir,
    name,h
    )

f=open("makeagravray.sh","w")
f.write("""#PBS -S /bin/bash
#PBS -j oe
#PBS -o log/%s.log
#PBS -t 1-%d
cd $PBS_O_WORKDIR
python makeagravray.py %s
"""%(runstr,nprocs,options))
f.close()

System("cp makeagravray.sh %s"%outdir)
