#!/usr/bin/env python
from gravray import *
from os import system
from sys import stderr,stdout

#############################################################
#USAGE
#############################################################
usage="""Make a GRT analysis of a whole geopgraphic region.

python makeagravray.py <date> <deg|rad> <file.geopos> <file.locals> <sname> <qvel> [<altitude>]

Where:

   <date>: date of analysis.
           Use format MM/DD/CCYY HH:MM:SS.dcm UTC-L.

   <deg|rad>: indicate if positions in geo.files are in deg or rad

   <file.geopos>: file with the geographical position (generated with
                  generatedirections.py)

   <file.locals>: file with the initial impacting conditions at each
                  site (generate using generatelocals.py)

   <sname>: name of the simulation

   <qvel>: use apex dependent velocities?

   <altitude>: altitude where the rays start in meters (default:
               80,000 m)

Output:

   Analysis directory : data/grt-CCYYMMDDHHMMSS-<id>

      where <id>: unique identifier depending on analysis configuration

   Location elements, ray-lat_<latitude>__lon_<longitude>.data: for
   each location it creates a file with the asymptotic orbit elements of the
   rays thrown from that site.

   Location probabilities,
   ray-lat_<latitude>__lon_<longitude>.data.prob: probability
   associated to each ray.

   geographic.prob: probability for each location. 

"""

#############################################################
#INPUT
#############################################################
try:
    iarg=1
    date=argv[iarg];iarg+=1
    deg=argv[iarg];iarg+=1
    geofile=argv[iarg];iarg+=1
    inifile=argv[iarg];iarg+=1
    QVEL=int(argv[iarg]);iarg+=1
    NAME=argv[iarg];iarg+=1
except:
    print usage
    exit(1)

#ALTITUDE
try:
    h=argv[iarg];iarg+=1
except:h=8e4

#############################################################
#PARAMETERS
#############################################################
#SET 1 IF YOU WANT THAT VELOCITY BE DEPENDENT OF QAPEX

#NAME OF THE ANALYSIS
#QVEL=0;NAME="Chalyabinsk apex dependent velocities (N=540 locals)"
#QVEL=1;NAME="Chalyabinsk average velocities (N=540 locals)"
#QVEL=0;NAME="Hawaii average velocities (N=540 locals)"
#QVEL=1;NAME="Hawaii apex dependent velocities (N=540 locals)"
#QVEL=0;NAME="Madagascar average velocities (N=540 locals)"
#QVEL=1;NAME="Madagascar apex dependent velocities (N=540 locals)"

#############################################################
#PROCESS INPUTS
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
makestr="QVEL=%d & NAME=%s"%(QVEL,NAME)
md5str=MD5STR(makestr,len=6)

#OUTPUT DIRECTORY
outdir="data/grt-%s-%s"%(tstring,md5str)
System("mkdir -p %s"%outdir)
f=open(outdir+"/.config","w")
f.write("%s\n%s\n%s\n"%(md5str,tstring,makestr))
System("cp %s %s/geographic.dat"%(geofile,outdir))
System("cp %s %s/locals.dat"%(inifile,outdir))
f.close()
print>>stderr,"Running analysis at %s..."%outdir

#LAT,LON UNITS
if deg=='deg':UNIT=1
else:UNIT=RAD

#DETERMINE NUMBER OF QAPEX
datos=np.loadtxt(inifile)
Ninitial=len(datos)
napex=datos.shape[1]-2

#############################################################
#GENERATE OBSERVERS MATRICES (ECLIPTIC AND APEX POSITION)
#############################################################
cmd="make && ./whereonearth.exe '%s'"%date
system(cmd)
system("cp -r scratch/observers-matrices.dat %s/"%outdir)

#############################################################
#MAKE ANALYSIS
#############################################################
#GEOGRAPHIC POSITIONS
data=np.loadtxt(geofile)
try:
    lons=data[:,2]*UNIT
    lats=data[:,3]*UNIT
    npoints=len(lats)
except:
    lons=[data[2]*UNIT]
    lats=[data[3]*UNIT]
    npoints=1

timeIt(stream=stderr)
print "Analysing %d points..."%npoints
f=open(outdir+"/geographic.prob","w")
for i in xrange(npoints):

    lat=lats[i]
    lon=lons[i]

    #==================================================
    #THROW RAYS FROM THIS PLACE
    #==================================================
    print>>stderr,"*"*80,"\nCalculating elements for location %d/%d: lat = %e, lon = %e...\n"%(i,npoints-1,lat,lon),"*"*80
    
    outfile="rays-lat_%.5e__lon_%.5e.data"%(lat,lon)
    print t,lat,lon,h,inifile,QVEL,outdir,outfile
    cmd="./throwrays.exe %.9e %.5e %.5e %.4e %s %d %s/%s"%(t,lat,lon,h,inifile,QVEL,outdir,outfile)
    print "Executing:",cmd
    system(cmd)
    timeIt(stream=stderr)

    #==================================================
    #CALCULATE PROBABILITIES
    #==================================================
    print "Calculating probabilities for this site"
    system("python analyseatsource.py %s/locals.dat %s/%s"%(outdir,outdir,outfile))
    timeIt(stream=stderr)

    #==================================================
    #CALCULATE TOTAL PROBABILITY
    #==================================================
    data=np.loadtxt("%s/%s.prob"%(outdir,outfile))
    p=data[:,7]
    Ptot=p.sum()/(1.0*Ninitial)
    print>>stderr,"Total probability:",Ptot
    f.write("%-+15.6f%-+15.6f%-+15.6f\n"%(lat,lon,Ptot))

f.close()
