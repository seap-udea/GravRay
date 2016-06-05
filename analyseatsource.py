from gravray import *

#############################################################
#INPUTS
#############################################################
fname=argv[1]
print "*"*80,"\nAnalysing data in '%s'\n"%fname,"*"*80

#############################################################
#CONSTANTS AND NUMERICAL PARAMETERS
#############################################################
Ninitial=530

#Weighting factor for computing density
def wFunction(d,h):
    """
    Schoenber B-spline function
    See: https://arxiv.org/pdf/1012.1885.pdf

    Plot:
    h=0.1
    sigma=wNormalization(h)
    fig=plt.figure()
    ax=fig.gca()
    ds=np.linspace(0,5*h,100)
    ws=np.array([sigma*wFunction(d,h) for d in ds])
    ax.plot(ds,ws)
    fig.savefig("scratch/weighting-shoenberg.png")

    Test it:
    from scipy.integrate import quad
    wnorm=lambda d:wFunction(d,h)*sigma
    print quad(wnorm,0,2*h)
    """
    q=d/h
    if q<1:w=0.25*(2-q)**3-(1-q)**3
    elif q<2:w=0.25*(2-q)**3
    else:w=0
    return w

def wNormalization(h):
    from scipy.integrate import quad
    sigma=1/quad(wFunction,0,2*h,args=(h,))[0]
    return sigma

#############################################################
#GET DATA FROM FILE
#############################################################
data=np.loadtxt(fname)
Norbits=data.shape[0]

Ncoll=Ninitial-Norbits

print "Number of initial conditions:",Ninitial
print "Number of succesfull orbits:",Norbits
print "Number of collisions:",Ncoll

#############################################################
#READ ELEMENTS OF IMPACTORS
#############################################################
qes=data[:,6]
ees=data[:,7]
ies=data[:,8]

cond=(ees<1)*(ies<180)

qes=qes[cond]
ees=ees[cond]
ies=ies[cond]

Nphys=ees.shape[0]

print "NUmber of bound, prograde orbits:",Nphys

aes=qes/(1-ees)

#############################################################
#NUMERICAL PARAMTERES
#############################################################

#Maximum weighted euclidean distance in configuration space
dmax=0.1

#Weighting function normalization
sigma=wNormalization(dmax)

#Normalization of number density
normal=1000.0*Ninitial

#############################################################
#COMPUTE DENSITY
#############################################################
Ptot=0

fp=open(fname+".prob","w")
for n in xrange(Nphys):
 
    q=qes[n]
    e=ees[n]
    i=ies[n]
    
    if (n%100)==0:
        print "Direction %d:"%n,q,e,i
        
    distform="POW(Perihelion_dist-%.17e,2)+POW(e-%.17e,2)+POW((i-%.17e)/45.,2)"%(q,e,i)
    distance=lambda qt,et,it:((qt-q)**2+(e-et)**2+((i-it)/90)**2)**0.5

    result=np.array(mysqlSelect("%s, Perihelion_dist, e, i"%distform,
                                "NEOS",
                                "where %s<%e order by %s"%(distform,(2*dmax)**2,distform),"array"))

    ntarg=result.shape[0]
    d2t,qt,et,it=0,0,0,0
    density=0
    if ntarg>0:
        d2t,qt,et,it=result[0,:]
        for target in result:
            d2,q,e,i=target
            d=d2**0.5
            p=sigma*wFunction(d,dmax)
            #print "\tTarget: d = %.2e, p = %.2e"%(d,p)
            density+=p

    Pn=density/normal
    Ptot+=Pn
    #print "Probability contribution: ",Pn/normal
    #print Pn
    fp.write("%+.3e %+.3e %+.3e %6d %.3e %.3e %.3e %+.5e\n"%(q,e,i,ntarg,qt,et,it,Pn))
    if n>100e2:break

fp.close()    
print "Total probability for this site: ",Ptot
