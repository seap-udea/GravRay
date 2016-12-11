#!/usr/bin/env python
from gravray import *
from os import path,system
from time import sleep

#############################################################
#USAGE
#############################################################
usage="""Make multiple Parallel GRT analysis of a whole geographic region.

   python multipleparallelgravray.py <nprocs> <dates_file> <deg|rad> <locations_file> <deg|rad> <directions_file> <velocities_file> [<height> <rundir>]

Where:

   <nprocs>: number of processors used for the analysis.

   <dates_file>: File the date of the analysis

   <deg|rad>: indicate if directions are in degrees or radians.

   <locations_file>: file with the geographical locatins (generated
                     using eg. generatedirections.py)

   <directions_file>: file with the impacting directons (generated
                     using eg. generatedirections.py)

   <velocities_file>: file with initial velocities (generated
                      using eg. generatevelocities.py)

   <altitude>: altitude where the rays start in meters (default:
               80,000 m)

   <rundir>: Directory with the launch scripts will be created.
"""

#############################################################
#PREPARE RUN
#############################################################
iarg=1
try:
    nprocs=int(argv[iarg]);iarg+=1

    datesfile=argv[iarg];iarg+=1

    degloc=argv[iarg];iarg+=1
    locfile=argv[iarg];iarg+=1

    degdir=argv[iarg];iarg+=1
    dirfile=argv[iarg];iarg+=1

    qvel=0
    velfile=argv[iarg];iarg+=1
except:
    print usage
    exit(1)

#ALTITUDE
try:h=float(argv[iarg]);iarg+=1
except:h=8e4

#DIRECTORY
try:rundir=argv[iarg];iarg+=1
except:rundir="."

#############################################################
#PREPARE
#############################################################
System("make")

#############################################################
#READ DATES FILE
#############################################################
f=open(datesfile)
dates=[]
for date in f:
    dates+=[date.replace("\n","")]

#############################################################
#CREATE RUNS
#############################################################
i=1
for date in dates:
    makefile="%s/makeagravray-%d.sh"%(rundir,i)
    name='Run %s'%date
    print "Creating parallel gravray %d in %s..."%(i,makefile)
    cmd="python parallelgravray.py %d '%s'  %s %s   %s %s   %s %s   '%s'   %f   %s   0"%(nprocs,date,degloc,locfile,degdir,dirfile,qvel,velfile,name,h,makefile)
    print "\tCommand: ",cmd
    system(cmd)
    i+=1

print

#NUMBER OF PROCESSES PER JOB
num=int(System("grep 'PBS -t' %s/makeagravray-1.sh |awk -F'-' '{print $3}'"%rundir))+1
print "Number of jobs per process: %d"%num
timeIt(stream=stderr)

#############################################################
#SUBMIT RUNS
#############################################################
print 

print "Running jobs sucessively..."
i=1
sleeptime=20
bunmd5=System("md5sum %s |awk '{print $1}'"%datesfile)
bunfile="data/bundle-%s.tar"%bunmd5
print "\tBundle file '%s'..."%bunfile

if path.isfile(bunfile):system("rm bunfile")
for date in dates:
    grtid='grt-'+System("grep 'PBS -o' %s/makeagravray-%d.sh |awk '{print $3}' |awk -F'log/' '{print $2}' |awk -F'.' '{print $1}'"%(rundir,i))
    print "Jobs for date %d '%s' grtid = '%s'..."%(i,date,grtid)
    out=System("qsub %s/makeagravray-%d.sh"%(rundir,i))
    jobid=out.split("[")[0]
    print "\tRunning job %s..."%jobid
    exect=0
    while True:
        cmd="grep '%s\[[0-9]*\]' /var/spool/torque/server_logs/* |grep 'cput=' |wc -l"%jobid
        ncompleted=int(System(cmd))
        print "\t\tJobs ncompleted after %d secs: %d"%(exect,ncompleted)
        if ncompleted==num:
            #cmd="python mapatsource.py data/grt-%s-%s %d %f %f"%(date,grtid,qmatrix,qlat,qlon)
            print "\t\tMapping probabilities..."
            cmd="python mapatsource.py data/%s"%(grtid)
            system(cmd)
            system("cp data/%s/probability-map-contour.png data/%s/probability-map-contour-%04d.png"%(grtid,grtid,i))
            print "\t\tPacking results..."
            System("rm data/%s/rays-*.data"%(grtid))
            system("tar rf %s -C data %s"%(bunfile,grtid))
            timeIt(stream=stderr)
            break
        else:
            #ESTIMATING REMAINING TIME AND ADJUSTING SLEEPTIME
            if ncompleted>0:sleept=sleeptime/(1+ncompleted/2.)
            else:sleept=sleeptime
            print "\t\tWaiting %d secs..."%sleeptime
            exect+=sleept
            sleep(sleept)
    print "\tJob terminated"
    i+=1

#############################################################
#COMPRESS RESULTS
#############################################################
system("gzip %s"%bunfile)
print "Multiple run completed."
timeIt(stream=stderr)
