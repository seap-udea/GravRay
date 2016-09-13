from gravray import *

#############################################################
#USAGE
#############################################################
usage="""Make a GRT analysis of a whole geopgraphic region.

python analyseatsource.py <file.locals> <file.elements>

Where:

   <file.locals>: file with initial conditions (azimuth, elevation,
                  velocities).

   <file.elements>: file with resulting elements after analysis with
                    throwrays.exe

Output:

   <file.elements>.prob: file containing the probability associated to
                         each ray.  Columns:

      #1:q       2:e        3:i        4:ntarg  5:qclose  6:eclose  7:iclose  8:probability
      +7.935e-01 +3.545e-01 +1.039e+01    634   6.277e-01 2.484e-01 6.888e+00 +9.10427e-03
   
   where ntarg is the number of objects in the database with values of
   the orbital elements close to that of the test particle;
   qclose,eclose,iclose are the elements of the closest object in the
   database to the test particle; probability is the "normalized"
   probability for this point.

Example:
   

"""

#############################################################
#INPUTS
#############################################################
try:
    iarg=1
    inifile=argv[iarg];iarg+=1
    elements=argv[iarg];iarg+=1
except:
    print usage
    exit(1)

print "*"*80,"\nAnalyzing data in '%s'\n"%elements,"*"*80

#############################################################
#CONSTANTS AND NUMERICAL PARAMETERS
#############################################################

#############################################################
#GET DATA FROM FILE
#############################################################
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

Nhyp=len(ees[ees>=1])
Nret=len(ies[ies>=180])
cond=(ees<1)*(ies<180)

qes=qes[cond]
ees=ees[cond]
ies=ies[cond]

Nphys=ees.shape[0]

print "Filter:"
print TAB,"Number of hyperbolic orbits:",Nhyp
print TAB,"Number of retrograde orbits:",Nret
print TAB,"Number of bound, prograde orbits:",Nphys

aes=qes/(1-ees)

#############################################################
#NUMERICAL PARAMTERES
#############################################################
#Debugging
verb=0
adv=0

#Maximum weighted euclidean distance in configuration space
dmax=0.1

#Weighting function normalization
sigma=wNormalization(dmax)

#Normalization of number density
normal=5000.0*Ninitial

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
    if verb:print "Test particle:",q,e,i

    if (n%(Nphys/10))==0 and adv:
        print "Direction %d:"%n,q,e,i
        
    distform=drummondDistance(q,e,i)

    result=np.array(mysqlSelect("%s, Perihelion_dist, e, i"%distform,
                                "NEOS",
                                "where %s<%e order by %s"%(distform,(2*dmax)**2,distform),"array"))

    ntarg=result.shape[0]

    if verb:print TAB,"Number of targets:",ntarg
    d2c,qc,ec,ic=0,0,0,0
    density=0
    if ntarg>0:
        d2c,qc,ec,ic=result[0,:]
        n=0
        for target in result:
            d2,q,e,i=target
            if verb:print TAB,"Target %d:"%n,q,e,i,", Distance:",d2
            d=d2**0.5
            p=sigma*wFunction(d,dmax)
            if verb:print TAB,"Target: d = %.2e, p = %.6e"%(d,p)
            if verb:raw_input()
            density+=p
            n+=1

    Pn=density/normal
    Ptot+=Pn
    if verb:print TAB,"Probability contribution: ",Pn/normal
    fp.write("%+.3e %+.3e %+.3e %6d %.3e %.3e %.3e %+.5e\n"%(q,e,i,ntarg,qc,ec,ic,Pn))
    if verb:raw_input()
    if n>100e2:break

fp.close()    
print "Total probability for this site: ",Ptot
timeIt()
