#!/usr/bin/env python
from gravray import *
from os import system

date=argv[1]
parts=date.split(" ")
dparts=parts[0].split("/")
tparts=parts[1].split(":")
tstring="%04d%02d%02d%02d%02d%02d"%(int(dparts[2]),int(dparts[0]),int(dparts[1]),
                                    int(tparts[0]),int(tparts[1]),int(tparts[2]))
edir="ensamble-%s"%tstring
out=System("./whattimeisit.exe '%s' ET > /dev/null"%date)
t=float(out.split("\n")[4])

points="scratch/points-r1.00e+01-filtered.data"
ensdir="data/%s"%edir
System("mkdir -p %s"%ensdir)

data=np.loadtxt(points)
lons=data[:,2]*180/np.pi
lats=data[:,3]*180/np.pi
npoints=len(lats)

"""
#15 February 2013, 3:20 UTC (Chelyabinsk event)
#t=4.141704340e+08

#12 hours after Chelyabinsk event
#t=4.142136340e+08

#6 months after Chelyabinsk event
t=4.298088000e+08

#30 June 1908, 2:30 UTC (Tunguska event)
#t=-2.887703160e+09
"""

h=8.234765e+04

f=open("%s/ensamble.dat"%ensdir,"w")

for i in xrange(npoints):

    lat=lats[i]
    lon=lons[i]

    print "*"*80,"\nCalculating elements for lat = %e, lon = %e...\n"%(lat,lon),"*"*80
    f.write("%e %e\n"%(lat,lon))
    cmd="time ./launchmany.exe %.9e %.5e %.5e %.4e initial.dat %s/elements-lat_%.5e__lon_%.5e.data"%(t,lat,lon,h,ensdir,lat,lon)
    system(cmd)

f.close()
