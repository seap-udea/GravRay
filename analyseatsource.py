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
