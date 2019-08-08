#!/usr/bin/env python
from gravray import *
from os import system,path
from sys import stderr,stdout
import string
import random

#############################################################
#USAGE
#############################################################
usage="""Make a GRT analysis of a whole geographic region.

   python makeahyperbolic.py <date> <deg|rad> <locations_file> <deg|rad> <directions_file> <velocities_file> <sname> <altitude> [<runid> <qrepeat>]

Where:

   <date>: date of analysis.
           Use format MM/DD/CCYY HH:MM:SS.dcm UTC-L.

   <deg|rad>: indicate if directions are in degrees or radians.

   <locations_file>: file with the geographical locatins (generated
                     using eg. generatedirections.py)

   <directions_file>: file with the impacting directons (generated
                     using eg. generatedirections.py)

   <velocities_file>: file with initial velocities (generated
                      using eg. generatevelocities.py)

   <sname>: name of the simulation

   <altitude>: altitude where the rays start in meters (default:
               80,000 m)

   <runid>: Run identifier.

   <deltat>: Backwards integration time.

Output:

   data/grt-CCYYMMDDHHMMSS-<id>: Analysis directory, where <id>:
   unique identifier depending on analysis configuration

   ray-lat_<latitude>__lon_<longitude>.data: for each location it
   creates a file with the asymptotic orbit elements of the rays
   thrown from that site.

   ray-lat_<latitude>__lon_<longitude>.data.prob: location
   probabilities, probability associated to all rays at a given
   location.

   geographic.prob: probability for each location. 

"""

#############################################################
#INPUT
#############################################################
iarg=1
try:
    date=argv[iarg];iarg+=1

    degloc=argv[iarg];iarg+=1
    locfile=argv[iarg];iarg+=1

    degdir=argv[iarg];iarg+=1
    dirfile=argv[iarg];iarg+=1

    qvel=0
    velfile=argv[iarg];iarg+=1

    name=argv[iarg];iarg+=1

    h=float(argv[iarg]);iarg+=1

    deltat=float(argv[iarg]);iarg+=1
except:
    print usage
    exit(1)

#RUN ID
try:
    runid="%03d"%int(argv[iarg]);iarg+=1
except:runid=""

#ALTITUDE
try:qrepeat=int(argv[iarg]);iarg+=1
except:qrepeat=0

#COMMAND
cmd=" ".join(argv)

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
makestr="cmd = %s & qvel=%d & name=%s & dirmd5=%s & velmd5=%s"%(cmd,qvel,name,dirmd5,velmd5)

md5str=MD5STR(makestr,len=6)

#OUTPUT DIRECTORY
outdir="data/grt-%s-%s-%s"%(name,tstring,md5str)
System("mkdir -p %s"%outdir)
f=open(outdir+"/.config","w")
f.write("%s\n%s\n%s\n"%(md5str,tstring,makestr))
System("cp %s %s/locations.dat"%(locfile,outdir))
System("cp %s %s/directions.dat"%(dirfile,outdir))
System("cp %s %s/velocities.dat"%(velfile,outdir))
f.close()
print>>stderr,"Running analysis at %s..."%outdir

#UNITS
UNITS=dict(deg=1,rad=RAD)
UNITLOC=UNITS[degloc]
UNITRAD=UNITS[degdir]

#DETERMINE NUMBER OF QAPEX
datos=np.loadtxt(velfile)
try:datos[:,0]
except:datos=np.array([datos])
napex=datos.shape[1]-1

#############################################################
#GENERATE INITIAL CONDITIONS
#############################################################
inifile="%s/initials%s.dat"%(outdir,runid)
f=open(inifile,"w")
f.write("#0:h\t1:A\t2:v\t3:qapex\n")
vels=np.loadtxt(velfile)
try:vels[:,0]
except:vels=np.array([vels])

datadir=np.loadtxt(dirfile)
try:
    bs=datadir[:,3]
    dirs=datadir[:,2:]
except:
    bs=datadir[3]
    dirs=np.array([datadir[2:]])
    
ndirs=datadir[bs>0].shape[0]

n=0
Azs=[]
Els=[]
nvel=len(vels[:,0])

for direction in dirs:
    l=np.mod(direction[0]*UNITRAD,360)
    b=direction[1]*UNITRAD
    if b<0:continue
    Azs+=[l]
    Els+=[b]
    for i in xrange(nvel):
        f.write("%-+20.4e%-+20.4e"%(b,l))
        for j in xrange(napex):
            f.write("%-+20.4e"%(vels[i,j]))
        f.write("%-+20d\n"%(-1))
        n+=1
f.close()
system("cp %s %s/initials.dat"%(inifile,outdir))
Ninitial=n
print "%d rays prepared..."%n

#############################################################
#MAKE ANALYSIS
#############################################################
data=np.loadtxt(locfile)
try:
    lons=data[:,2]*UNITLOC
    lats=data[:,3]*UNITLOC
    npoints=len(lats)
except:
    lons=[data[2]*UNITLOC]
    lats=[data[3]*UNITLOC]
    npoints=1

timeIt(stream=stderr)
print "Analysing %d locations..."%npoints

ranstr=''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))

for i in xrange(npoints):

    lat=lats[i]
    lon=lons[i]

    #==================================================
    #THROW RAYS FROM THIS PLACE
    #==================================================
    print>>stderr,"*"*80,"\nCalculating elements for location %d/%d: lat = %e, lon = %e...\n"%(i+1,npoints,lat,lon),"*"*80
    
    outfile="rays-lat_%.5e__lon_%.5e.data"%(lat,lon)
    if not path.isfile("%s/%s"%(outdir,outfile)) or qrepeat:
        qvel=0
        cmd="./throwlongrays.exe %d %.9e %.5e %.5e %.4e %s %d %s/%s %e"%(n,t,lat,lon,h,inifile,
                                                                         qvel,outdir,outfile,deltat)
        print "Executing: %s"%cmd
        system(cmd)
