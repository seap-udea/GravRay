#!/usr/bin/env python
from gravray import *
from os import system
from sys import stderr,stdout

#############################################################
#USAGE
#############################################################
usage="""Make a GRT analysis of a whole geographic region.

   python makeagravray.py <date> <deg|rad> <locations_file> <deg|rad> <directions_file> <qvel> <velocities_file> <sname> [<altitude>]

Where:

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
    inifile=argv[iarg];iarg+=1

    qvel=int(argv[iarg]);iarg+=1
    velfile=argv[iarg];iarg+=1

    name=argv[iarg];iarg+=1
except:
    print usage
    exit(1)

#ALTITUDE
try:
    h=argv[iarg];iarg+=1
except:h=8e4

#############################################################
#INPUTS ANALYSIS
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
makestr="qvel=%d & name=%s"%(qvel,name)
md5str=MD5STR(makestr,len=6)

#OUTPUT DIRECTORY
outdir="data/grt-%s-%s"%(tstring,md5str)
System("mkdir -p %s"%outdir)
f=open(outdir+"/.config","w")
f.write("%s\n%s\n%s\n"%(md5str,tstring,makestr))
System("cp %s %s/locations.dat"%(locfile,outdir))
System("cp %s %s/directions.dat"%(inifile,outdir))
System("cp %s %s/velocities.dat"%(velfile,outdir))
f.close()
print>>stderr,"Running analysis at %s..."%outdir

#UNITS
UNITS=dict(deg=1,rad=RAD)
UNITLOC=UNITS[degloc]
UNITRAD=UNITS[degdir]

#DETERMINE NUMBER OF QAPEX
datos=np.loadtxt(inifile)
Ninitial=len(datos)
napex=datos.shape[1]-2

#############################################################
#GENERATE OBSERVERS MATRICES (ECLIPTIC AND APEX POSITION)
#############################################################
system("./whereonearth.exe '%s'"%date)
system("cp -r scratch/observers-matrices.dat %s/"%outdir)

exit(0)

#############################################################
#MAKE ANALYSIS
#############################################################
#GEOGRAPHIC POSITIONS
data=np.loadtxt(geofile)
try:
    lons=data[:,2]*UNITLOC
    lats=data[:,3]*UNITLOC
    npoints=len(lats)
except:
    lons=[data[2]*UNITLOC]
    lats=[data[3]*UNITLOC]
    npoints=1

timeIt(stream=stderr)
print "Analysing %d points..."%npoints

import string
import random
ranstr=''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
f=open(outdir+"/geographic-%s.prob"%ranstr,"w")
for i in xrange(npoints):

    lat=lats[i]
    lon=lons[i]

    #==================================================
    #THROW RAYS FROM THIS PLACE
    #==================================================
    print>>stderr,"*"*80,"\nCalculating elements for location %d/%d: lat = %e, lon = %e...\n"%(i,npoints-1,lat,lon),"*"*80
    
    outfile="rays-lat_%.5e__lon_%.5e.data"%(lat,lon)
    cmd="./throwrays.exe %.9e %.5e %.5e %.4e %s %d %s/%s"%(t,lat,lon,h,inifile,qvel,outdir,outfile)
    print "Executing:",cmd
    system(cmd)
    timeIt(stream=stderr)

    #==================================================
    #CALCULATE PROBABILITIES
    #==================================================
    print "Calculating probabilities for this site"
    cmd="python analyseatsource.py %s/locals.dat %s/%s"%(outdir,outdir,outfile)
    print cmd
    system(cmd)
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
