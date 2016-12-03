from gravray import *
from os import path,system
from sys import stderr,stdout

#############################################################
#USAGE
#############################################################
usage="""Update impact probabilities using a different model

python analyseatsource.py <grt.dir>

Where:

   <grt.dir>: directory with the complete grt analysis.

Output:

   Ray files with probabilities <rays>.data.prob

   geographic.prob

Example:
   
  python updateatsource.py data/grt-20130215032034-3CAA5C

"""

#############################################################
#INPUTS
#############################################################
try:
    iarg=1
    grtdir=argv[iarg];iarg+=1
except:
    print usage
    exit(1)

if not path.isdir(grtdir):
    print "Directory '%s' does not exist."%grtdir
    print 
    print usage
    exit(1)

print "*"*80,"\nUpdating probabilities in '%s'\n"%grtdir,"*"*80

#############################################################
#GET ALL RAY FILES
#############################################################
rays=System("ls -m %s/rays-*.data"%grtdir).split(",\n")
Nrays=len(rays)
inifile="%s/locals.dat"%grtdir
data=np.loadtxt(inifile)
Ninitial=len(data)

system("mv %s/geographic*.prob %s/geographic.prob.save"%(grtdir,grtdir))

i=1
f=open(grtdir+"/geographic.prob","w")
for ray in rays:

    print "#"*70,"\nUpdating %d/%d location\n"%(i,Nrays),"#"*70

    #GET LATITUDE AND LONGITUDE
    pcoords=ray.split("-lat")[1].split("__")
    lat=float(pcoords[0].replace("_",""))
    lon=float(pcoords[1].replace("lon_","").replace(".data",""))

    #RE ANALYZE PROBABILITY
    cmd="python analyseatsource.py %s %s"%(inifile,ray)
    print "Executing:",cmd
    system(cmd)

    #COMPUTE PROBABILITY FOR THIS DIRECTION
    data=np.loadtxt("%s.prob"%ray)
    p=data[:,7]
    Ptot=p.sum()/(1.0*Ninitial)
    
    #SHOW TOTAL
    print>>stderr,"Total probability:",Ptot
    f.write("%-+15.6f%-+15.6f%-+15.6f\n"%(lat,lon,Ptot))

    i+=1
f.close()
