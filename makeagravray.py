#!/usr/bin/env python
from gravray import *
from os import system

#############################################################
#USAGE
#############################################################
usage="""Make a GRT analysis of a whole geopgraphic region.

python makeagravray.py <date> <file.geopos> <file.initials> [<altitude>]

Where:

   <date>: date of analysis.  Use format MM/DD/CCYY HH:MM:SS.dcm
   UTC-L.

   <file.geopos>: file with the geographical position (generated with
   generatedirections.py)

   <file.initials>: file with the initial impacting conditions at each
   site (generate using generateinitialconditions.py)

   <altitude>: altitude where the rays start in meters (default:
   80,000 m)

"""

#############################################################
#INPUT
#############################################################
try:date=argv[1]
except:
    print usage
    exit(1)

#DATE ANALYSIS
parts=date.split(" ")
dparts=parts[0].split("/")
tparts=parts[1].split(":")
tstring="%04d%02d%02d%02d%02d%02d"%(int(dparts[2]),int(dparts[0]),int(dparts[1]),
                                    int(tparts[0]),int(tparts[1]),int(tparts[2]))
out=System("./whattimeisit.exe '%s' ET > /dev/null"%date)
t=float(out.split("\n")[4])

#FILES
try:
    geofile=argv[2]
    inifile=argv[3]
except:
    print usage
    exit(1)

#ALTITUDE
try:h=argv[4]
except:h=8e4

#OUTPUT DIRECTORY
outdir="data/grt-%s"%tstring
System("mkdir -p %s"%outdir)

#############################################################
#MAKE ANALYSIS
#############################################################
#GEOGRAPHIC POSITIONS
data=np.loadtxt(geofile)
lons=data[:,2]*180/np.pi
lats=data[:,3]*180/np.pi
npoints=len(lats)

timeIt()
print "Analysing %d points..."%npoints
for i in xrange(npoints):

    lat=lats[i]
    lon=lons[i]

    #THROW RAYS FROM THIS PLACE
    print "*"*80,"\nCalculating elements for lat = %e, lon = %e...\n"%(lat,lon),"*"*80
    cmd="./throwrays.exe %.9e %.5e %.5e %.4e %s %s/rays-lat_%.5e__lon_%.5e.data"%(t,lat,lon,h,inifile,outdir,lat,lon)
    system(cmd)
    timeIt()

    #ANALYSE PROBABILITY OF THESE RAYS

    break
