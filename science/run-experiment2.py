from gravray import *
from os import system

sims=[
    'QVEL=0;NAME="Chelyabinsk average velocities (N=540 locals)";site="chelyabinsk"',
    'QVEL=0;NAME="Madagascar average velocities (N=540 locals)";site="madagascar"',
    'QVEL=0;NAME="Hawaii average velocities (N=540 locals)";site="hawaii"',
    #'QVEL=1;NAME="Chelyabinsk apex dependent velocities (N=540 locals)";site="chelyabinsk"',
    #'QVEL=1;NAME="Hawaii apex dependent velocities (N=540 locals)";site="hawaii"',
    #'QVEL=1;NAME="Madagascar apex dependent velocities (N=540 locals)";site="madagascar"',
    #'QVEL=0;NAME="Antartica average velocities (N=540 locals)";site="antartica"',
    #'QVEL=1;NAME="Antartica apex dependente velocities (N=540 locals)";site="antartica"',
    #'QVEL=1;NAME="Australia apex dependente velocities (N=540 locals)";site="australia"',
    #'QVEL=1;NAME="Atlantic apex dependente velocities (N=540 locals)";site="atlantic"',
    ]

#MAKE STRING
for sim in sims:
    exec(sim)
    makestr="QVEL=%d & NAME=%s"%(QVEL,NAME)
    md5str=MD5STR(makestr,len=6)
    odir="data/grt-20130215032034-%s"%md5str
    system("mkdir -p %s"%odir)
    cmd="python makeagravray.py '02/15/2013 03:20:34 UTC' deg geographic.dat.%s locals.dat %d '%s' > %s/grt.log"%(site,QVEL,NAME,odir)
    print "Running:",cmd
    system(cmd)
