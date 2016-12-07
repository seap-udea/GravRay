#!/usr/bin/env python
from gravray import *
from os import system

#############################################################
#USAGE
#############################################################
usage="""Make a GRT analysis of a whole geographic region.

   python backupgrt.py <depth> <grt_directory>

Where:

   <depth>: 1 for only probability data, 2 for all data.

"""

#############################################################
#INPUT
#############################################################
iarg=1
depth=argv[iarg];iarg+=1

try:
    depth=int(depth)
except:
    print usage
    exit(1)

#############################################################
#COMPILE
#############################################################
system("rm scratch/results*.*")

outdirs=",".join(argv[2:])
print "Packing directories %s at depth %d..."%(outdirs,depth)

for outdir in argv[2:]:
    print "\tPacking %s..."%outdir
    if depth==1:
        system("tar rf scratch/results.tar %s/{probability.prob,initials.dat,observers-matrices.dat}"%(outdir))
    else:
        system("tar rf scratch/results.tar %s/"%(outdir))

system("gzip scratch/results.tar")
