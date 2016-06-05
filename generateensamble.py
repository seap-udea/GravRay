#!/usr/bin/env python
import numpy as np
from sys import exit
from os import system

points="scratch/points-r1.00e+01-filtered.data"

data=np.loadtxt(points)
lons=data[:,2]*180/np.pi
lats=data[:,3]*180/np.pi
npoints=len(lats)

t=4.141704340e+08
h=8.234765e+04

f=open("data/ensamble/ensamble.dat","w")
for i in xrange(npoints):

    lat=lats[i]
    lon=lons[i]

    print "*"*80,"\nCalculating elements for lat = %e, lon = %e...\n"%(lat,lon),"*"*80
    f.write("%e %e\n"%(lat,lon))
    cmd="time ./launchmany.exe %.9e %.5e %.5e %.4e initial.dat data/ensamble/elements-lat_%.5e__lon_%.5e.data"%(t,lat,lon,h,lat,lon)
    system(cmd)

f.close()
