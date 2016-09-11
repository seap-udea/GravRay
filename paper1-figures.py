from gravray import *
from scipy import signal as sig
from scipy.optimize import curve_fit

#############################################################
#INPUTS
#############################################################
def histOutline(histIn,binsIn):
    stepSize = binsIn[1] - binsIn[0]

    bins = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    data = np.zeros(len(binsIn)*2 + 2, dtype=np.float)
    for bb in range(len(binsIn)):
        bins[2*bb + 1] = binsIn[bb]
        bins[2*bb + 2] = binsIn[bb] + stepSize
        if bb < len(histIn):
            data[2*bb + 1] = histIn[bb]
            data[2*bb + 2] = histIn[bb]

    bins[0] = bins[1]
    bins[-1] = bins[-2]
    data[0] = 0
    data[-1] = 0

    return (bins, data)

def cumDistrib(x,h,x0=None,xn=None):
    norm=h.sum()
    F=np.zeros_like(h)
    F[0]=h[0]
    for i in xrange(1,h.shape[0]):
        F[i]=F[i-1]+h[i]
    F=F/(1.0*norm)
    if x0 is not None:
        x=np.concatenate(([x0],x))
        F=np.concatenate(([0],F))
    if xn is not None:
        x=np.concatenate((x,[xn]))
        F=np.concatenate((F,[1]))
    Fcum=np.vstack((x,F)).transpose()
    return Fcum

VESC=11.1
def theoVimp(v,p,a,vso):
    pv=p*((v-VESC)/vso)*np.exp(-((v-VESC)/vso)**a)
    return pv

#############################################################
#CONSTANTS AND NUMERICAL PARAMETERS
#############################################################
REARTH=6371 #km

FIGDIR="paper1-figures/"

#############################################################
#READ DATA
#############################################################
qload=0
try:qload=int(argv[1])
except:pass

if qload>0:
    props="Perihelion_dist, Aphelion_dist, e, i"
    condition="where Perihelion_dist<=1 and Perihelion_dist>=(1-e)/(1+e) and (e>0 and e<1)"
    listdict=mysqlSelect(props,
                         condition=condition)
    elements=listdict2matrix(listdict,keys=props.split(", "))

#############################################################
#FUNCTIONS
#############################################################
def Ncum(E):
    """
    Number of events per year with energy E (kton) or larger
    """
    #Brown et al. (2013)
    a=3.31
    b=-0.68
    N=a*E**b
    return N

def pExponential(t,lamb):
    """
    Probability density that time between two random events happen be
    t given that the rate of events be lamb
    """
    p=lamb*np.exp(-lamb*t)
    return p

def cumpExponential(T,lamb):
    """
    Probability that the time between two events be less than t given a
    rate of events lamb
    """
    P=1-np.exp(-lamb*T)
    return P
    

def numpExponential(N,T,lamb):
    """
    Probability that N events (with rate lamb) happen in time T
    """
    from scipy.misc import factorial
    P=(lamb*T)**N*np.exp(-lamb*T)/factorial(N)
    return P

#############################################################
#ROUTINES
#############################################################
def distanceTunguskaChelyabisnk():
    #Tunguska coordinates
    pT=s2d(101,57)*DEG,s2d(60,55)*DEG

    #Chelyabinsk
    pC=55.150*DEG,64.410*DEG

    #Angular distance
    qdist=arcDistance(pT,pC)
    print "Angular distance: %.2f deg (%.5f rad)"%(qdist*RAD,qdist)

    #Linear distance
    dist=REARTH*qdist
    print "Linear distance: %.2f km"%dist

    #Solid angle
    sangle=2*PI*(1-np.cos(qdist))
    print "Solid angle: %.5f pi"%(sangle/PI)

    #Fraction of the sphere
    fsphere=sangle/(4*PI)
    print "Fraction of the sphere: %.6f"%fsphere

    #Rate of impacts
    lamb2013=Ncum(500)
    Dt=1/lamb2013
    print "Rate larger than Chelyabinsk (Brown+2013):",lamb2013
    print "Average time between impactes:",Dt
    
    #Rate according Harris (2012)
    lamb2012=7.3e-3
    print "Rate larger than Chelyabinsk (Harris+2013):",lamb2012

    #Probability that time between impacts be less than 20 years
    P=cumpExponential(20.0,lamb2012)
    print "Probability < 20 years (Harris+2012):",P
    P=cumpExponential(20.0,lamb2013)
    print "Probability < 20 years (Brown+2013):",P

    #Probability that time between impacts be 20 years
    p=pExponential(20.0,lamb2012)
    print "Probability = 20 years:",p

    #Probability that 1 event happen in 20 years
    p1=numpExponential(1,20,lamb2012)
    print "Probability than 1 event happen in 20 years (Harris+2012):",p1

    #Probability that 1 to 30 events happen in time T
    Ns=np.arange(30)
    print Ns

def plotGeographicPositions():
    data=np.loadtxt("directions-r5.00e+00.dat")
    ssgood=data[:,2:]
    N=len(ssgood)

    plt.close("all")

    print "Map of Poisson Noise points..."
    fig=plt.figure(figsize=(8,8))

    print "Map of Random Noise points..."
    #ax=fig.add_axes([0.05,0.6,0.9,0.5])
    ax=fig.add_subplot(212)
    
    #Phong-like distribution 
    #See: http://www.cs.rutgers.edu/~decarlo/readings/mcrt-sg03c.pdf, p. 23
    #For a uniform distribution over a sphere use n=0
    n=0
    r1=rand(N)
    r2=rand(N)
    q=np.arccos((1-2*r1)**(1/(n+1.)))*RAD
    f=2*np.pi*r2*RAD
    l=90.-q

    proj='moll'
    map=drawMap(proj=proj,proj_opts=dict(lon_0=180,ax=ax),
                pars=[-45,0,45],mers=[0,45,90,135,180,225,270,315],
                pars_opts=dict(labels=[1,1,0,0],fontsize=8),
                mers_opts=dict(labels=[0,0,0,1],fontsize=8))
    plotMap(map,f,l,lw=0,
            marker='o',color='r',ms=4,mec='none')
    ax.set_title(r"Random sampling: $N = %d$"%N,position=(0.5,1.05))

    #ax=fig.add_axes([0.05,0.05,0.9,0.5])
    ax=fig.add_subplot(211)
    proj='moll'
    radius=5.0
    map=drawMap(proj=proj,proj_opts=dict(lon_0=180,ax=ax),
                pars=[-45,0,45],mers=[0,45,90,135,180,225,270,315],
                pars_opts=dict(labels=[1,1,0,0],fontsize=8),
                mers_opts=dict(labels=[0,0,0,1],fontsize=8))
    plotMap(map,np.mod(ssgood[:,0]*RAD,360),ssgood[:,1]*RAD,lw=0,
            marker='o',color='r',ms=4,mec='none')
    ax.set_title(r"Poisson sampling: $N = %d,\;r = 5^\circ$"%N,position=(0.5,1.05))


    fig.tight_layout()
    fig.savefig(FIGDIR+"InitialGeographicConditions.png")

def apexVelocityDistribution():

    """
    Run after veclocityDistribution
    """
    plt.close("all")

    # ########################################
    # PARAMETERS
    # ########################################
    ifreq=1
    perc=90
    pnorm=1

    # ########################################
    # READ OBJECTS DATA
    # ########################################
    data=np.loadtxt(FIGDIR+"objects.data")
    
    # IMPACT INCOMING DIRECTION
    tapexs=data[:,21]
    ntapex=len(tapexs)
    tmin=tapexs.min()
    tmax=tapexs.max()

    print "Number of simulated objects: ",ntapex
    
    # ########################################
    # RANGE OF INPUT DIRECTIONS
    # ########################################
    tetas=[(0,60),(120,180),(60,120)]
    ntetas=len(tetas)
    
    # ########################################
    # AVERAGE SPEED DISTRIBUTION
    # ########################################
    print "Speed distribution for the whole sample"
    nbins=50
    vimps=data[:,17]

    print "\tPercentile %.2f:"%perc,np.percentile(vimps,perc)
    print "\tWidth:",2*vimps.std()

    # HISTOGRAM
    h,v=np.histogram(vimps,nbins,(VESC,4*VESC))
    vb=(v[:-1]+v[1:])/2

    # NORMALIZATION
    if pnorm:h=h/(1.0*h.max())
    else:h=h/(1.0*ntapex)
    bins,n=histOutline(h,vb)
    label="Average"

    # HISTOGRAM FIT
    p=0.1;a=0.1;vso=1.0
    fit=curve_fit(theoVimp,vb,h,(p,a,vso))
    param=fit[0]
    print "\tFitting parameters: ",param

    # ANALYTICAL FREQUENCY
    vbs=np.linspace(VESC,vb.max(),100)
    pvbs=theoVimp(vbs,*param)

    # ANALYTICAL CUMULATIVE
    Fcum=cumDistrib(vbs,pvbs)
    np.savetxt("vdistrib-qapex_0_180.data",Fcum)

    # NORMALIZATION
    if pnorm:ppl=pvbs/pvbs.max()
    else:ppl=pvbs

    # ########################################
    # PLOT
    # ########################################
    fig=plt.figure()
    ax=fig.gca()
    
    ax.plot(vbs,ppl,label=label,color='k',lw=10,alpha=0.2,zorder=-10)

    # ########################################
    # DISTRIBUTION CALCULATED PER INTERVAL
    # ########################################

    # ====================
    # PARAMETERS
    # ====================
    colors=['b','r','g']
    ncolors=3
    nbins=50

    for i in xrange(ntetas):

        # ====================
        # RANGE OF THETAS
        # ====================
        tm,tM=tetas[i]
        cond=(tapexs>=tm)*(tapexs<tM)
        nvs=len(tapexs[cond])

        print "Interval %d: tm,tM,n="%i,tm,tM,nvs

        # ====================
        # SPEEDS
        # ====================
        vimps=data[cond,17]
        label=r"$%.0f^{\circ}\leq\theta_{\rm apex}<%.0f^\circ$"%(tm,tM)

        print "\tPercentile %.0f:"%perc,np.percentile(vimps,perc)
        print "\tWidth:",2*vimps.std()

        # ====================
        # HISTOGRAM
        # ====================
        h,v=np.histogram(vimps,nbins,(VESC,4*VESC))
        vb=(v[:-1]+v[1:])/2

        # NORMALIZATION
        if pnorm:h=h/(1.0*h.max())
        else:h=h/(1.0*nvs)

        #PLOT 
        bins,n=histOutline(h,vb)
        #ax.plot(bins,n,color=colors[i%ncolors],lw=2)

        # ====================
        # FIT
        # ====================
        try:
            p=0.1;a=0.1;vso=1.0
            fit=curve_fit(theoVimp,vb,h,(p,a,vso))
            param=fit[0]
            print "\tFitting parameteres: ",param

            # PLOTTING
            vbs=np.linspace(VESC,vb.max(),100)
            pvbs=theoVimp(vbs,*param)

            # ANALYTICAL CUMULATIVE
            Fcum=cumDistrib(vbs,pvbs)
            np.savetxt("vdistrib-qapex_%.0f_%.0f.data"%(tm,tM),Fcum)

            if pnorm:ppl=pvbs/pvbs.max()
            else:ppl=pvbs
            ax.plot(vbs,ppl,label=label,color=colors[i%ncolors],lw=2)

        except:
            print "\tFit did not converge..."

           
    # ########################################
    # DECORATION
    # ########################################
    ax.set_xlim((11,40))
    ax.set_yticks([])
    ax.set_xlabel("$v$ (km/s)",fontsize=14)
    ax.set_ylabel("Scaled frequency",fontsize=14)
    ax.legend(loc='upper right')

    fig.savefig(FIGDIR+"ApexVelocityDistributions.png")

    # ########################################
    # PLOT DISTRIBUTION OF DIRECTIONS
    # ########################################
    print "Distribution of directions..."
    fig=plt.figure()
    ax=fig.gca()
    ax.hist(tapexs,40)
    fig.savefig(FIGDIR+"DirectionDistributions.png")
    
    # ########################################
    # DELETE DATA
    # ########################################
    del(data)

def velocityDistribution(el,verbose=0):

    from scipy.optimize import curve_fit

    nbodies=el.shape[0]

    #Atmospheric entry height
    H=80.0 #km

    #OBRITAL VELOCITY
    vorb=np.sqrt(6.67E-11*2E30/1.5E11)/1E3
    print "Earth orbital velocity (km/s): ",vorb

    #ESCAPE VELOCITY
    vesc2=(2*5.98E24/2E30)*(1.5e8/(6.471e3+H))
    vesc=np.sqrt(vesc2)*vorb
    print "Earth escape velocity (v^2,v,v[km/s]): ",vesc2,np.sqrt(vesc2),vesc

    #SUN ESCAPE VELOCITY
    vescsun=np.sqrt((2*6.67e-11*2E30)/1.5e11)/1e3
    print "Sun escape velocity (v[km/s]): ",vescsun

    allindexes=np.arange(nbodies)
    nsamples=10
    nsubsample=int(nbodies*0.8)
    nbins=40

    hmeanimps=np.zeros(nbins)
    hmeaninfinis=np.zeros(nbins)
    hmeaninfsrcs=np.zeros(nbins)

    fo=open(FIGDIR+"objects.data","w")
    fo.write("%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s%-25s\n"%(\
                                         '#1:q','2:Q','3:a','4:p','5:e','6:i','7:O','8:w','9:wpf',
                                         '10:xdot','11:ydot','12:zdot','13:va*vorb',
                                         '14:xdotrel','15:ydotrel','16:zdotrel','17:vinf*vorb',
                                         '18:vimp*vorb',
                                         '19:xapex','20:yapex','21:zapex','22:tapex'))
    np.random.seed(4)
    for n in xrange(nsamples):

        indexes=np.random.choice(allindexes,nsubsample)

        velimps=[]
        velinfinis=[]
        velinfsrcs=[]

        for i in indexes:

            #ORBITAL ELEMENTS
            q=el[i,0]
            Q=el[i,1]
            e=el[i,2]
            i=el[i,3]

            #DERIVED PROPERTIES
            a=q/(1-e)
            p=q*(1+e)
            hop=1/np.sqrt(p);

            #ASCENDING OR DESCENDING NODE
            O=0
            if np.random.rand()>0.5:O=180*DEG
            
            #PERIHELIA ARGUMENT
            if O==0:
                w=np.arccos((p-1)/e)
                wpf=0.0*DEG
            else:
                w=np.arccos((1-p)/e)
                wpf=180*DEG

            if verbose:
                print "Mock asteroid:"
                print "\t","q,Q,a,p,e,i = ",q,Q,a,p,e,i
                print "\t","Omega = ",O*RAD
                print "\t","omega = ",w*RAD

            #ASTEROID VELOCITY (VIS VIVA)
            r=1.0
            va2=(2/r-1/a)
            va=np.sqrt(va2)

            #VELOCITY COMPONENTS
            xdot=-hop*e*np.cos(O)*np.sin(w)
            ydot=+hop*np.cos(O)*np.cos(i*DEG)*(np.cos(wpf)+e*np.cos(w))
            zdot=+hop*np.sin(i*DEG)*(np.cos(wpf)+e*np.cos(w))
            vadot2=xdot**2+ydot**2+zdot**2
            
            if verbose:
                print "3-D Velocity:"
                print "\t","xdot,ydot,zdot = ",xdot,ydot,zdot
                print "\t","Speed (vis viva) = ",va2
                print "\t","Speed (vector) = ",vadot2

            #RELATIVE VELOCITY: vrel = vearth - vast
            xdotrel=xdot-0
            ydotrel=ydot-1
            zdotrel=zdot-0
            vinf2=xdotrel**2+ydotrel**2+zdotrel**2
            vinf=np.sqrt(vinf2)

            if verbose:
                print "3-D Relative velocity:"
                print "\t","xdotrel,ydotrel,zdotrel = ",xdotrel,ydotrel,zdotrel
                print "\t","Relative speed (velocity at infinite) = ",vinf

            #IMPACT VELOCITY
            vimp2=vinf2+vesc2
            vimp=np.sqrt(vimp2)
            
            if verbose:
                print "Impact velocity:"
                print "\t","vimp^2 = ",vimp2
                print "\t","vimp = ",vimp
                print "\t","vimp [km/s] = ",vimp*vorb

            #APEX DIRECTION
            xapex=+xdotrel
            yapex=+zdotrel
            zapex=-ydotrel
            tapex=np.arccos(zapex/vinf)
            fapex=np.arctan2(yapex,xapex)

            if verbose:
                print "Apex direction:"
                print "\t","xapex,yapex,zapex = ",xapex,yapex,zapex
                print "\t","Theta apex = ",tapex*RAD

            fo.write("%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e%-+25.17e\n"%(
                                         q,Q,a,p,e,i,O,w,wpf,
                                         xdot,ydot,zdot,va*vorb,
                                         xdotrel,ydotrel,zdotrel,vinf*vorb,
                                         vimp*vorb,
                                         xapex,yapex,zapex,tapex*RAD))
            if verbose:
                raw_input()
                
            # STORE VALUES
            velinfsrcs+=[va*vorb]
            velinfinis+=[vinf*vorb]
            velimps+=[vimp*vorb]

        velimps=np.array(velimps)
        velinfsrcs=np.array(velinfsrcs)
        velinfinis=np.array(velinfinis)

        h,ximp=np.histogram(velimps,nbins,(vesc,4*vesc))
        hmeanimps+=h

        h,xinfini=np.histogram(velinfinis,nbins,(0,4*vesc))
        hmeaninfinis+=h

        h,xinfsrc=np.histogram(velinfsrcs,nbins,(0.5*vorb,vescsun))
        hmeaninfsrcs+=h

    fo.close()
    hmeanimps/=1.0*nsamples
    hmeaninfinis/=1.0*nsamples
    hmeaninfsrcs/=1.0*nsamples

    ###################################################
    #MAIN PLOT
    ###################################################

    fig=plt.figure()
    ax=fig.gca()

    """
    bins,n=histOutline(hmeaninfsrcs/nsubsample,xinfsrc)
    ax.plot(bins,n,label=r'Source $v_\infty$')
    """

    bins,n=histOutline(hmeaninfinis/nsubsample,xinfini)
    ax.plot(bins,n,label='Initial velocities $v_\infty$')

    bins,n=histOutline(hmeanimps/nsubsample,ximp)
    ax.plot(bins,n,label=r'Impact velocities $v_{\rm imp}$')
    
    #THEORETICAL FIT
    """
    def theoVimp(v,p,a,vso):
        pv=p*(vs-vesc)*np.exp(-(vs-vesc)**a/vso)
        return pv

    vs=np.concatenate(([vesc],(ximp[:-1]+ximp[1:])/2))
    pvs=np.concatenate(([0],hmeanimps/nsubsample))

    p=0.1;a=0.1;vso=1.0
    fit=curve_fit(theoVimp,(vs,),pvs,(p,a,vso))
    """
    def theoVimp(v,p,a,vso):
        pv=p*(v-vesc)*np.exp(-(v-vesc)**a/vso)
        return pv
    vs=np.concatenate(([vesc],(ximp[:-1]+ximp[1:])/2))
    pvs=np.concatenate(([0],hmeanimps/nsubsample))

    p=0.1;a=0.1;vso=1.0
    fit=curve_fit(theoVimp,vs,pvs,(p,a,vso))
    param=fit[0]
    print param[2]**(1/param[1])
    print "Fit result:",param
    ax.plot(vs,theoVimp(vs,*param),'k',label=r'Continuous fit: $\alpha=0.88,\;\Delta v=3.04$ km/s')

    #MEAN VALUE
    vmean=((ximp[:-1]+ximp[1:])/2*hmeanimps/nsubsample).sum()
    ax.axvline(vmean,color='k',ls='dashed')
    ax.text(vmean+1,0.005,r"$<v_{\rm imp}>=%.1f$ km/s"%vmean)
    
    #==============================
    #DECORATION
    #==============================
    ax.axvspan(0,vesc,color='k',alpha=0.2)
    ax.set_xlabel(r"$v$ (km/s)",fontsize=14)
    ax.set_ylabel(r"Frequency",fontsize=14)
    ax.legend(loc="upper right",fontsize=10)
    ax.set_xlim((0,42))
    ax.set_ylim((0,0.09))

    fig.savefig(FIGDIR+"VelocityDistribution.png")

    ###################################################
    #SAVE CUMMULATIVE IMPACT VELOCITY DISTRIBUTION
    ###################################################
    x=ximp
    xbin=(ximp[:-1]+ximp[1:])/2
    hmean=hmeanimps

    Fmean=np.zeros_like(hmean)
    for i in xrange(1,hmean.shape[0]):
        Fmean[i]=Fmean[i-1]+hmean[i-1]
    Fmean/=nsubsample
    np.savetxt(FIGDIR+"velocity-cum.data",np.vstack(
        (np.append(xbin,x[-1]),
         np.append(Fmean,[1]))
    ).transpose())

    ###################################################
    #SAVE HISTOGRAM
    ###################################################
    np.savetxt(FIGDIR+"velocity-pdf.data",np.vstack((xbin,hmean/(1.0*nsubsample))).transpose())

def readECSS2008():
    data=np.loadtxt("dev/ECSS.data")
    x=np.concatenate(data.transpose()[::2,:])
    n=np.concatenate(data.transpose()[1::2,:])

#############################################################
#EXECUTE
#############################################################
#distanceTunguskaChelyabisnk()
#plotGeographicPositions()
