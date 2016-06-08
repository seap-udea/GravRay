#!/usr/bin/env python
import numpy as np
from sys import exit
from os import system

points="scratch/points-r1.00e+01-filtered.data"
ensdir="data/ensamble3"

data=np.loadtxt(points)
lons=data[:,2]*180/np.pi
lats=data[:,3]*180/np.pi
npoints=len(lats)

#15 February 2013, 3:20 UTC (Chelyabinsk event)
#t=4.141704340e+08

#12 hours after Chelyabinsk event
#t=4.142136340e+08

#6 months after Chelyabinsk event
t=4.298088000e+08

#30 June 1908, 2:30 UTC (Tunguska event)
#t=-2.887703160e+09

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
