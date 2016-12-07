#!/usr/bin/env python
from gravray import *
from os import system

#############################################################
#USAGE
#############################################################
usage="""Make a GRT analysis of a whole geographic region.

   python backupgrt.py <grt_directory>

"""

#############################################################
#INPUT
#############################################################
iarg=1
try:
    outdirs=argv[iarg];iarg+=1
except:
    print usage
    exit(1)

#############################################################
#COMPILE
#############################################################
system("rm scratch/results*.*")

outdirs=",".join(argv[1:])
print "Packing directories %s..."%outdirs

for outdir in argv[1:]:
    print "\tPacking %s..."%outdir
    system("tar rf scratch/results.tar %s/{probability.prob,initials.dat,observers-matrices.dat}"%(outdir))

system("gzip scratch/results.tar")
