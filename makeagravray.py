#!/usr/bin/env python
from gravray import *
from os import system

#############################################################
#USAGE
#############################################################
usage="""Make a GRT analysis of a whole geopgraphic region.

python makeagravray.py <date> <deg|rad> <file.geopos> <file.locals> [<altitude>]

Where:

   <date>: date of analysis.
           Use format MM/DD/CCYY HH:MM:SS.dcm UTC-L.

   <deg|rad>: indicate if positions in geo.files are in deg or rad

   <file.geopos>: file with the geographical position (generated with
                  generatedirections.py)

   <file.locals>: file with the initial impacting conditions at each
                  site (generate using generatelocals.py)

   <altitude>: altitude where the rays start in meters (default:
               80,000 m)

Output:

   

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
QVEL=1

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

#OUTPUT DIRECTORY
outdir="data/grt-%s"%tstring
System("mkdir -p %s"%outdir)

#LAT,LON UNITS
if deg=='deg':UNIT=1
else:UNIT=RAD

#DETERMINE NUMBER OF QAPEX
datos=np.loadtxt(inifile)
napex=datos.shape[1]-2

#############################################################
#MAKE ANALYSIS
#############################################################
#GEOGRAPHIC POSITIONS
data=np.loadtxt(geofile)
lons=data[:,2]*UNIT
lats=data[:,3]*UNIT
npoints=len(lats)

timeIt()
print "Analysing %d points..."%npoints
for i in xrange(npoints):

    lat=lats[i]
    lon=lons[i]

    #THROW RAYS FROM THIS PLACE
    print "*"*80,"\nCalculating elements for lat = %e, lon = %e...\n"%(lat,lon),"*"*80
    cmd="./throwrays.exe %.9e %.5e %.5e %.4e %s %d %s/rays-lat_%.5e__lon_%.5e.data"%(t,lat,lon,h,inifile,QVEL,outdir,lat,lon)
    print "Executing:",cmd
    system(cmd)
    timeIt()

    #ANALYSE PROBABILITY OF THESE RAYS

    break
