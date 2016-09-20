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
inifile="%s/locals.dat"%grtdir
Ninitial=len(inifile)

f=open(grtdir+"/geographic.prob.test","w")
for ray in rays:

    #GET LATITUDE AND LONGITUDE
    
    exit(0)

    #RE ANALYZE PROBABILITY
    cmd="python analyseatsource.py %s %s"%(inifile,ray)
    system(cmd)

    #COMPUTE PROBABILITY FOR THIS DIRECTION
    data=np.loadtxt("%s.prob"%ray)
    p=data[:,7]
    Ptot=p.sum()/(1.0*Ninitial)
    
    #SHOW TOTAL
    print>>stderr,"Total probability:",Ptot
    f.write("%-+15.6f%-+15.6f%-+15.6f\n"%(lat,lon,Ptot))
    break

exit(0)

initials=np.loadtxt(inifile)
Ninitial=len(initials)

data=np.loadtxt(elements)
Norbits=data.shape[0]
Ncoll=Ninitial-Norbits

print "Basic properties:"
print TAB,"Number of initial conditions:",Ninitial
print TAB,"Number of succesfull orbits:",Norbits
print TAB,"Number of collisions:",Ncoll

#############################################################
#READ ELEMENTS OF IMPACTORS
#############################################################
qes=data[:,9]
ees=data[:,10]
ies=data[:,11]
qxs=data[:,15]
aes=qes/(1-ees)

#NEW-INCLUDING NODE AND PERIHELION
Omegas=data[:,12]
omegas=data[:,13]

Nhyp=len(ees[ees>=1])
Nret=len(ies[ies>=180])
cond=(ees<1)*(ies<180)*(aes<40)

qes=qes[cond]
ees=ees[cond]
ies=ies[cond]
aes=aes[cond]
qxs=qxs[cond]

#NEW-INCLUDING NODE AND PERIHELION
Omegas=Omegas[cond]
omegas=omegas[cond]

Nphys=ees.shape[0]

print "Filter:"
print TAB,"Number of hyperbolic orbits:",Nhyp
print TAB,"Number of retrograde orbits:",Nret
print TAB,"Number of bound, prograde orbits:",Nphys

#############################################################
#NUMERICAL PARAMTERES
#############################################################
#Debugging
verb=0
adv=0

#Maximum weighted euclidean distance in configuration space
#dmax=0.1
#NEW
dmax=0.15

#Weighting function normalization
sigma=wNormalization(dmax)

#Maximim value of the smoothing kernel
wmax=sigma*wFunction(0,dmax)

#Normalization of number density
#normal=2000.0
normal=1.0 #NEW 

#Flux function parameters
#Obtained with paper1-figures, apexVelocityDistribution()
fparam=(0.9721768,6.84870896,2.40674371)

#############################################################
#COMPUTE DENSITY
#############################################################
Ptot=0

timeIt()
fp=open(elements+".prob","w")
for n in xrange(Nphys):
 
    q=qes[n]
    e=ees[n]
    i=ies[n]

    #NEW
    Omega=Omegas[n]
    omega=omegas[n]
    a=aes[n]

    qx=qxs[n]
    flux=theoFlux_DoubleTrigCos(qx,*fparam)

    if verb:print "Test particle:",q,e,i

    if (n%(Nphys/10))==0 and adv:
        print "Direction %d:"%n,q,e,i
        
    #distform=drummondDistance(q,e,i)
    #NEW
    distform=zappalaDistance(a,e,np.sin(i*DEG),Omega,omega)

    """
    result=np.array(mysqlSelect("%s, Perihelion_dist, e, i"%distform,
                                "NEOS",
                                "where %s<%e order by %s"%(distform,(2*dmax)**2,distform),"array"))
    """
    #NEW
    result=np.array(mysqlSelect("%s, Perihelion_dist, e, i, sini, a, Node, Peri"%distform,
                                "NEOS",
                                "where %s<%e order by %s"%(distform,(2*dmax)**2,distform),"array"))

    ntarg=result.shape[0]
    """
    if ntarg>1000:
        print "Elements:",a,e,np.sin(i*DEG),Omega,omega
        raw_input()
    """

    if verb:print TAB,"Number of targets:",ntarg
    
    d2,qt,et,it,sinit,at,Ot,ot=0,0,0,0,0,0,0,0

    density=0
    if ntarg>0:
        n=0

        #NEW
        d2,qt,et,it,sinit,at,Ot,ot=result[0,:]

        for target in result:
            #d2,q,e,i=target
            #NEW
            d2,qt,et,it,sinit,at,Ot,ot=target
            d=d2**0.5
            p=sigma*wFunction(d,dmax)
            if verb:print "q=%.3f,%.3f"%(q,qt),"e=%.3f,%.3f"%(e,et),"i=%.3f,%.3f"%(i,it),"sini=%.3f,%.3f"%(np.sin(i*DEG),sinit),"a=%.3f,%.3f"%(a,at),"O = %.3f,%.3f"%(Omega,Ot),"o = %.3f,%.3f"%(omega,ot),"d = %.3f"%d,"p=",p
            density+=p
            n+=1
        if verb:print "Density:",density
    else:
        density=0

    Pu=density/wmax
    Pn=flux*Pu

    if verb:print TAB,"Probability contribution: ",Pn
    Ptot+=Pn/normal
    fp.write("%+.3e %+.3e %+.3e %6d %.3e %.3e %.3e %+.5e %.5e %.2f %.3e %.3e %.3e\n"%(q,e,i,ntarg,qt,et,it,Pn/normal,Pu/normal,qx,at,Ot,ot))

    if verb:raw_input()
    if n>100e2:break

fp.close()    

#Normalize total probability
Ptot=Ptot/(1.0*Ninitial)
print "Total probability for this site: ",Ptot
timeIt()
