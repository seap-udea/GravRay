"""
This file let you to plan
"""
from gravray import *
from os import system
from sys import argv

nsim=int(argv[1])

if nsim==1:
    sim='QVEL=0;NAME="Whole World Tunguska (N=540 locals)";site="Whole world"'
    exec(sim)
    makestr="QVEL=%d & NAME=%s"%(QVEL,NAME)
    md5str=MD5STR(makestr,len=6)
    odir="data/grt-19080630001400-%s"%md5str
    system("mkdir -p %s"%odir)
    cmd="python makeagravray.py '06/30/1908 00:14:00 UTC' rad util/data/directions-r7.00e+00.data locals.dat %d '%s' > %s/grt.log"%(QVEL,NAME,odir)
    print "Running:",cmd
    system(cmd)

if nsim==2:
    sim='QVEL=0;NAME="Whole World before Chelyabinsk (N=540 locals)";site="Whole world"'
    exec(sim)
    makestr="QVEL=%d & NAME=%s"%(QVEL,NAME)
    md5str=MD5STR(makestr,len=6)
    odir="data/grt-20130214212034-%s"%md5str
    system("mkdir -p %s"%odir)
    cmd="python makeagravray.py '02/14/2013 21:20:34 UTC' rad util/data/directions-r7.00e+00.data locals.dat %d '%s' > %s/grt.log"%(QVEL,NAME,odir)
    print "Running:",cmd
    system(cmd)

if nsim==3:
    sim='QVEL=0;NAME="Whole World Africa Event (N=540 locals)";site="Whole world"'
    exec(sim)
    makestr="QVEL=%d & NAME=%s"%(QVEL,NAME)
    md5str=MD5STR(makestr,len=6)
    odir="data/grt-19630803164500-%s"%md5str
    system("mkdir -p %s"%odir)
    cmd="python makeagravray.py '08/03/1963 16:45:00 UTC' rad util/data/directions-r7.00e+00.data locals.dat %d '%s' > %s/grt.log"%(QVEL,NAME,odir)
    print "Running:",cmd
    system(cmd)

if nsim==4:
    sim='QVEL=0;NAME="Whole World Chelyabinsk (N=540 locals)";site="Whole world"'
    exec(sim)
    makestr="QVEL=%d & NAME=%s"%(QVEL,NAME)
    md5str=MD5STR(makestr,len=6)
    odir="data/grt-20130214212034-%s"%md5str
    system("mkdir -p %s"%odir)
    cmd="python makeagravray.py '02/15/2013 03:20:34 UTC' rad util/data/directions-r7.00e+00.data locals.dat %d '%s' > %s/grt.log"%(QVEL,NAME,odir)
    print "Running:",cmd
    system(cmd)

if nsim==5:
    sim='QVEL=1;NAME="Whole World Chelyabinsk (High resolution)";site="Whole world"'
    exec(sim)
    makestr="QVEL=%d & NAME=%s"%(QVEL,NAME)
    md5str=MD5STR(makestr,len=6)
    odir="data/grt-20130214212034-%s"%md5str
    system("mkdir -p %s"%odir)
    cmd="python makeagravray.py '02/15/2013 03:20:34 UTC' rad geographic.big.1 util/data/locals-r1.00e+01-v50.data %d '%s' > %s/grt.log"%(QVEL,NAME,odir)
    print "Running:",cmd
    system(cmd)

if nsim==6:
    sim='QVEL=1;NAME="Whole World Chelyabinsk (High resolution)";site="Whole world"'
    exec(sim)
    makestr="QVEL=%d & NAME=%s"%(QVEL,NAME)
    md5str=MD5STR(makestr,len=6)
    odir="data/grt-20130214212034-%s"%md5str
    system("mkdir -p %s"%odir)
    cmd="python makeagravray.py '02/15/2013 03:20:34 UTC' rad geographic.big.2 util/data/locals-r1.00e+01-v50.data %d '%s' > %s/grt.log"%(QVEL,NAME,odir)
    print "Running:",cmd
    system(cmd)

if nsim==7:
    sim='QVEL=1;NAME="Whole World Chelyabinsk (High resolution)";site="Whole world"'
    exec(sim)
    makestr="QVEL=%d & NAME=%s"%(QVEL,NAME)
    md5str=MD5STR(makestr,len=6)
    odir="data/grt-20130214212034-%s"%md5str
    system("mkdir -p %s"%odir)
    cmd="python makeagravray.py '02/15/2013 03:20:34 UTC' rad geographic.big.3 util/data/locals-r1.00e+01-v50.data %d '%s' > %s/grt.log"%(QVEL,NAME,odir)
    print "Running:",cmd
    system(cmd)

if nsim==8:
    sim='QVEL=1;NAME="Whole World Chelyabinsk (High resolution)";site="Whole world"'
    exec(sim)
    makestr="QVEL=%d & NAME=%s"%(QVEL,NAME)
    md5str=MD5STR(makestr,len=6)
    odir="data/grt-20130214212034-%s"%md5str
    system("mkdir -p %s"%odir)
    cmd="python makeagravray.py '02/15/2013 03:20:34 UTC' rad geographic.big.4 util/data/locals-r1.00e+01-v50.data %d '%s' > %s/grt.log"%(QVEL,NAME,odir)
    print "Running:",cmd
    system(cmd)

