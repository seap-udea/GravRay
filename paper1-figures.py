from gravray import *
from scipy import signal as sig
from scipy.optimize import curve_fit
from os import system,path
import csv
from spiceypy import wrapper as spy

#############################################################
#INPUTS
#############################################################
plt.close("all")

VESC=11.1
def theoVimp(v,p,a,vso):
    pv=p*((v-VESC)/vso)*np.exp(-((v-VESC)/vso)**a)
    return pv

def theoFlux_Trig(q,f,a,b):
    fv=f*(q*DEG)**a*np.sin(q*DEG)**b
    return fv

def theoFlux_DoubleTrigSin(q,fa,fb,a,b):
    if q<=90:
        fv=fa*np.sin(q*DEG)**a
    else:
        fv=fb*np.sin(q*DEG)**b
    return fv
theoFlux_DoubleTrigSin=np.vectorize(theoFlux_DoubleTrigSin)

def theoFlux_DoubleTrigCos(q,f,a,b):
    if q<=90:
        fv=f*np.cos(np.pi/2-q*DEG)**a
    else:
        fv=f*np.cos(np.pi/2-q*DEG)**b
    return fv
theoFlux_DoubleTrigCos=np.vectorize(theoFlux_DoubleTrigCos)

def theoFlux_DoubleExp(q,f,a,qa,b,qb):
    if qa<0 or qb<0:return 1e100
    if q<=90:
        fv=np.exp(-((90-q)*DEG/qa)**a)
    else:
        fv=((q-90)*DEG/qb)**(-b)
        #fv=np.exp(-((q-90)*DEG/qb)**b)
    fv*=f
    return fv
theoFlux_DoubleExp=np.vectorize(theoFlux_DoubleExp)

def theoFlux_PowExp(q,f,a,b,qt):
    fv=f*(q*DEG)**a/(1+np.exp((q**DEG/qt)**b))
    return fv

def softSignal(s,window=3):
    ns=len(s)
    ss=s
    for i in xrange(window,ns-window):
        ies=np.arange(i-window,i+window,1)
        ss[i]=s[ies].mean()
    return ss

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

if qload==1:
    print "Getting elements for closest NEAs..."
    props="Perihelion_dist, Aphelion_dist, e, i"
    condition="where Perihelion_dist<=1 and Perihelion_dist>=(1-e)/(1+e) and (e>0 and e<1)"
    listdict=mysqlSelect(props,
                         condition=condition)
    elements=listdict2matrix(listdict,keys=props.split(", "))
    print len(elements)," objects discovered..."

elif qload==2:
    print "Getting elements for all NEAs..."
    props="Perihelion_dist, e, i, sin(i*PI()/180), a, Node, Peri"
    condition="where NEO_flag and e>0 and e<1 and Perihelion_dist<=1"
    listdict=mysqlSelect(props,
                         condition=condition)
    elements=listdict2matrix(listdict,keys=props.split(", "))
    print len(elements)," objects discovered..."

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
    Fcum[-1,1]=1.0
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
            Fcum[-1,1]=1.0
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
    h,t=np.histogram(tapexs,40)
    tm=(t[1:]+t[:-1])/2

    #AREA CORRECTION
    hc=h/np.cos(np.pi/2-tm*DEG)

    #REGULARIZATION
    hc=softSignal(hc,window=4)
    hc=hc/hc.max()

    #FIT ANALYTICAL FUNCTION
    function=theoFlux_DoubleTrigCos
    f=0.1;a=1;b=2
    fit=curve_fit(function,tm,hc,(f,a,b))
    param=fit[0]
    print "Flux parameters:",param
    qts=np.linspace(0,180.0,100)
    fts=function(qts,*param)

    #PLOT
    fig=plt.figure()
    ax=fig.gca()
    bins,n=histOutline(hc,t)
    ax.plot(bins,n,label='Softed scaled flux')
    ax.plot(qts,fts,'k',lw=2,label='Trigonometric fit')
    ax.set_xlim((0,180))

    ax.set_xlabel(r"$\theta_{\rm apex}^\circ$",fontsize=14)
    ax.set_ylabel("Scaled frequency",fontsize=14)
    ax.legend(loc='upper left',fontsize=12)

    fig.savefig(FIGDIR+"DirectionDistributions.png")

    # ########################################
    # DELETE DATA
    # ########################################
    del(data)

def velocityDistribution(el,verbose=0):

    from scipy.optimize import curve_fit

    ###################################################
    #OBSERVED FIREBALLS
    ###################################################
    data=np.loadtxt("util/data/fireballs.data")
    fvimps=data[:,5]
    cond=fvimps!=123456789
    fvimps=fvimps[cond]
    ndataf=len(fvimps)
    hf,xf=np.histogram(fvimps,10)
    xmedf=(xf[:-1]+xf[1:])/2
    herrf=np.sqrt(hf)
    normf=1.0*ndataf*(xf[1]-xf[0])
    hnf=hf/normf
    herrnf=herrf/normf
    

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
    
    # FIREBALL DISTRIBUTION
    # binsf,nf=histOutline(hf,xf)
    # ax.plot(binsf,nf,'k-',lw=3,label='Observed fireballs velocities')
    ax.errorbar(xmedf,hnf,yerr=herrnf,fmt='o',label='NASA NEO program')
    

    # INITIAL DISTRIBUTION
    bins,n=histOutline(hmeaninfinis/nsubsample,xinfini)
    ax.plot(bins,n,label='Initial velocities $v_\infty$')

    # IMPACT DISTRIBUTION
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
    ax.plot(vs,theoVimp(vs,*param),'k',label=r'Continuous fit:$\alpha=0.88,\;\Delta v=3.04$ km/s')

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
    ax.legend(loc="upper right",fontsize=10,frameon=False)
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
    data=np.vstack((x,n)).transpose()
    np.savetxt("util/data/ECSS.data",data)

def plotKernel():
    h=0.1
    sigma=wNormalization(h)
    fig=plt.figure()
    ax=fig.gca()
    ds=np.linspace(0,5*h,100)
    ws=np.array([sigma*wFunction(d,h) for d in ds])
    ax.plot(ds,ws)
    figfile=FIGDIR+"WeightingShoenberg.png"
    print "Plotting ",figfile
    fig.savefig(figfile)

def testParticle():
    system("make && python throwaray.py 54.456093 63.492323 8e+04 17.664618 104.975030 -2.045864e+01 '02/15/2013 3:20:34 UTC' -2.0")

def showDistrib(el):

    #==================================================
    #DATA FROM NEA
    #==================================================    
    #LOAD DATA
    qes=el[:,0];qmin=qes.min();qmax=qes.max()
    ees=el[:,1];emin=ees.min();emax=ees.max()
    ies=np.log10(el[:,2]);imin=ies.min();imax=ies.max()
    aes=qes/(1-ees**2);amin=aes.min();amax=aes.max()

    #INTERVALS
    qlow=qmin;qup=qmax
    alow=amin;aup=amax*0+3
    elow=emin;eup=emax*0
    ilow=imin*0;iup=imax

    #==================================================
    #CALCULATE DENSITY
    #==================================================    
    bins=30

    print "Showing Source Distribution..."

    #==================================================
    #PLOT
    #==================================================    
    #cmap='gray'
    #interpolation='hanning'
    interpolation='nearest'
    cmap='rainbow'
    cmap='jet'

    factor=1.0
    pprop=dict(ms=5,mec='none')

    #fig=plt.figure(figsize=(18,7))
    fig=plt.figure()
    fsize=18

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # q vs. e
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    condi=np.abs(qes)>=0
    H,xe,ye=np.histogram2d(qes[condi],ees[condi],bins=bins,normed=True)
    axae=fig.add_subplot(131)
    scale=(qup-qlow)/(eup-elow)
    img=axae.imshow(H.transpose(),origin='lower',
                  interpolation=interpolation,
                  extent=(qmin,qmax,emin,emax),aspect=scale/factor,cmap=cmap)
    axae.set_xlabel("$q$ (AU)",fontsize=fsize)
    axae.set_ylabel("$e$",fontsize=fsize)
    axae.set_xlim((qlow,qup))
    axae.set_ylim((elow,eup))
    
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # q vs. i
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    condi=np.abs(qes)>=0
    H,xe,ye=np.histogram2d(qes[condi],ies[condi],bins=bins,normed=True)
    axai=fig.add_subplot(132)
    scale=(qup-qlow)/(iup-ilow)
    img=axai.imshow(H.transpose(),origin='lower',
                    interpolation=interpolation,
                    extent=(qmin,qmax,imin,imax),aspect=scale/factor,cmap=cmap)
    axai.set_xlabel("$q$ (AU)",fontsize=fsize)
    axai.set_ylabel("$\log(i^\circ)$",fontsize=fsize)
    axai.set_xlim((qlow,qup))
    axai.set_ylim((ilow,iup))
    axai.set_title("Marginal distributions of %d NEAs"%len(el),position=(0.5,1.05),fontsize=20)    

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # e vs. i
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    condi=np.abs(qes)>=0
    H,xe,ye=np.histogram2d(ees[condi],ies[condi],bins=bins,normed=True)
    axei=fig.add_subplot(133)
    scale=(eup-elow)/(iup-ilow)
    img=axei.imshow(H.transpose(),origin='lower',
                    interpolation=interpolation,
                    extent=(emin,emax,imin,imax),aspect=scale/factor,cmap=cmap)
    axei.set_xlabel("$e$",fontsize=fsize)
    axei.set_ylabel("$\log(i^\circ)$",fontsize=fsize)

    axei.set_xlim((elow,eup))
    axei.set_ylim((ilow,iup))

    fig.tight_layout()
    fig.savefig(FIGDIR+"NEODistribution.png")

def pointMap(el,fname,sname,title=None):
    """
    Run this with option 2 in command line
    """
    #==================================================
    #DATA FROM NEA
    #==================================================    
    #LOAD NEOS DATA
    qes=el[:,0];qmin=qes.min();qmax=qes.max()
    ees=el[:,1];emin=ees.min();emax=ees.max()
    ies=np.log10(el[:,2]);imin=ies.min();imax=ies.max()
    aes=qes/(1-ees**2);amin=aes.min();amax=aes.max()

    #==================================================
    #DATA FOR SITE
    #==================================================    
    data=np.loadtxt(fname)
    qdata=data[:,9]
    edata=data[:,10]
    idata=np.log10(data[:,11])
    cond=edata<1
    edata=edata[cond]
    qdata=qdata[cond]
    idata=idata[cond]
    adata=qdata/(1-edata)

    #==================================================
    #CHELYABINSK IMPACT
    #==================================================    
    eimp=0.6
    iimp=np.log10(6.0)
    qimp=0.75
    
    #INTERVALS
    qlow=qmin;qup=qmax
    elow=emin;eup=emax
    ilow=imin*0;iup=imax

    #==================================================
    #CALCULATE DENSITY
    #==================================================    
    bins=30

    print "Showing Source Distribution for site '%s'..."%sname

    #==================================================
    #PLOT
    #==================================================    
    interpolation='nearest'
    cmap='rainbow'
    #cmap='jet'

    factor=1.0
    pprop=dict(ms=8,mec='none')
    iprop=dict(ms=15,mec='k',color='w')

    fig=plt.figure(figsize=(18,7))
    fsize=18

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # q vs. e
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    condi=np.abs(qes)>=0
    H,xe,ye=np.histogram2d(qes[condi],ees[condi],bins=bins,normed=True)
    axae=fig.add_subplot(131)
    scale=(qup-qlow)/(eup-elow)
    img=axae.imshow(H.transpose(),origin='lower',
                  interpolation=interpolation,
                  extent=(qmin,qmax,emin,emax),aspect=scale/factor,cmap=cmap)
    axae.set_xlabel("$q$ (AU)",fontsize=fsize)
    axae.set_ylabel("$e$",fontsize=fsize)
    axae.plot(qdata,edata,'ko',**pprop)
    axae.plot([qimp],[eimp],'kv',**iprop)

    axae.set_xlim((qlow,qup))
    axae.set_ylim((elow,eup))

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # q vs. i
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    condi=np.abs(qes)>=0
    H,xe,ye=np.histogram2d(qes[condi],ies[condi],bins=bins,normed=True)
    axai=fig.add_subplot(132)
    scale=(qup-qlow)/(iup-ilow)
    img=axai.imshow(H.transpose(),origin='lower',
                    interpolation=interpolation,
                    extent=(qmin,qmax,imin,imax),aspect=scale/factor,cmap=cmap)
    axai.set_xlabel("$q$ (AU)",fontsize=fsize)
    axai.set_ylabel("$\log(i^\circ)$",fontsize=fsize)
    axai.plot(qdata,idata,'ko',**pprop)
    axai.plot([qimp],[iimp],'kv',**iprop)
    axai.set_xlim((qlow,qup))
    axai.set_ylim((ilow,iup))

    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # e vs. i
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    condi=np.abs(qes)>=0
    H,xe,ye=np.histogram2d(ees[condi],ies[condi],bins=bins,normed=True)
    axei=fig.add_subplot(133)
    scale=(eup-elow)/(iup-ilow)
    img=axei.imshow(H.transpose(),origin='lower',
                    interpolation=interpolation,
                    extent=(emin,emax,imin,imax),aspect=scale/factor,cmap=cmap)
    axei.set_xlabel("$e$",fontsize=fsize)
    axei.set_ylabel("$\log(i^\circ)$",fontsize=fsize)
    axei.plot(edata,idata,'ko',**pprop)    
    axei.plot([eimp],[iimp],'kv',**iprop)

    axei.set_xlim((elow,eup))
    axei.set_ylim((ilow,iup))

    if title is None:
        title=sname
    axai.text(0.5,1.05,"Point Distribution for site %s"%title,transform=axai.transAxes,fontsize=20,ha='center')    

    fig.tight_layout()
    fig.savefig(FIGDIR+"pointMap-%s.png"%sname)
    

def allPoints(el):
    """
    Generate point distributions:
    """
    #grtid="078645" #120 test particles
    grtid="9104FD" #271 test particles
    pointMap(el,'data/grt-20130215032034-%s/rays-lat_5.44000e+01__lon_6.35000e+01.data'%grtid,'Chelyabinsk')
    pointMap(el,'data/grt-20130215032034-%s/rays-lat_-1.89000e+01__lon_4.75000e+01.data'%grtid,'Madagascar')
    pointMap(el,'data/grt-20130215032034-%s/rays-lat_2.18000e+01__lon_-1.57000e+02.data'%grtid,'Hawaii')
    pointMap(el,'data/grt-20130215032034-%s/rays-lat_-7.78500e+01__lon_1.66400e+02.data'%grtid,'Antartica')

def analyzeProbability():

    #"""
    #HAWAII
    loc="lat_2.18000e+01__lon_-1.57000e+02"
    grtid="ED2863" #APEX
    date="20130215032034"
    cmd="python analyseatsource.py data/grt-%s-%s/locals.dat data/grt-%s-%s/rays-%s.data"%(date,grtid,date,grtid,loc)
    print "Executing:",cmd
    system(cmd)
    #"""

    #"""
    #MADAGASCAR
    loc="lat_-1.89000e+01__lon_4.75000e+01"
    grtid="21A94B" #APEX
    date="20130215032034"
    cmd="python analyseatsource.py data/grt-%s-%s/locals.dat data/grt-%s-%s/rays-%s.data"%(date,grtid,date,grtid,loc)
    print "Executing:",cmd
    system(cmd)
    #"""

    #"""
    #CHELYABINSK
    loc="lat_5.44000e+01__lon_6.35000e+01"
    grtid="FA0C86" #APEX 
    date="20130215032034"
    cmd="python analyseatsource.py data/grt-%s-%s/locals.dat data/grt-%s-%s/rays-%s.data"%(date,grtid,date,grtid,loc)
    print "Executing:",cmd
    system(cmd)
    #"""

def experimen2(el):

    """
    Run this with option 2 in command line
    """

    grtid="7F5C81" 
    loc="lat_-1.50000e+01__lon_1.35000e+02"
    print "Probability:",np.loadtxt("data/grt-20130215032034-%s/geographic.prob"%grtid)
    pointMap(el,'data/grt-20130215032034-%s/rays-%s.data'%(grtid,loc),'Australia_apex',title='Australia (apex speed distribution)')

    data=np.loadtxt("data/grt-20130215032034-%s/rays-%s.data"%(grtid,loc))
    qapex=np.unique(data[:,15])
    print qapex
    print np.percentile(qapex,[0,25,50,90])

    grtid="7368BB" 
    loc="lat_1.50000e+01__lon_-4.50000e+01"
    print "Probability:",np.loadtxt("data/grt-20130215032034-%s/geographic.prob"%grtid)
    pointMap(el,'data/grt-20130215032034-%s/rays-%s.data'%(grtid,loc),'Atlantic_apex',title='Atlantic (apex speed distribution)')

    data=np.loadtxt("data/grt-20130215032034-%s/rays-%s.data"%(grtid,loc))
    qapex=np.unique(data[:,15])
    print qapex
    print np.percentile(qapex,[0,25,50,90])

    """
    grtid="19C993" 
    print "Probability:",np.loadtxt("data/grt-20130215032034-%s/geographic.prob"%grtid)
    pointMap(el,'data/grt-20130215032034-%s/rays-lat_5.44000e+01__lon_6.35000e+01.data'%grtid,'Chelyabinsk_average',title='Chelyabinsk (average speed distribution)')
    grtid="FA0C86" 
    print "Probability:",np.loadtxt("data/grt-20130215032034-%s/geographic.prob"%grtid)
    pointMap(el,'data/grt-20130215032034-%s/rays-lat_5.44000e+01__lon_6.35000e+01.data'%grtid,'Chelyabinsk_apex',title='Chelyabinsk (apex dependent speed distribution)')
    #"""

    """
    loc="lat_2.18000e+01__lon_-1.57000e+02"
    grtid="2BAF32" #
    print "Probability:",np.loadtxt("data/grt-20130215032034-%s/geographic.prob"%grtid)
    pointMap(el,'data/grt-20130215032034-%s/rays-%s.data'%(grtid,loc),'Hawaii_average',title='Hawaii (average speed distribution)')
    grtid="ED2863" #
    print "Probability:",np.loadtxt("data/grt-20130215032034-%s/geographic.prob"%grtid)
    pointMap(el,'data/grt-20130215032034-%s/rays-%s.data'%(grtid,loc),'Hawaii_apex',title='Hawaii (apex dependent speed distribution)')
    #"""

    """
    loc="lat_-1.89000e+01__lon_4.75000e+01"
    grtid="0F13D5" #
    print "Probability:",np.loadtxt("data/grt-20130215032034-%s/geographic.prob"%grtid)
    pointMap(el,'data/grt-20130215032034-%s/rays-%s.data'%(grtid,loc),'Madagascar_average',title='Madagascar (average speed distribution)')
    grtid="21A94B" #
    print "Probability:",np.loadtxt("data/grt-20130215032034-%s/geographic.prob"%grtid)
    pointMap(el,'data/grt-20130215032034-%s/rays-%s.data'%(grtid,loc),'Madagascar_apex',title='Madagascar (apex dependent speed distribution)')
    #"""

    """
    loc="lat_-7.78500e+01__lon_1.66400e+02"
    grtid="1C8476" #
    print "Probability:",np.loadtxt("data/grt-20130215032034-%s/geographic.prob"%grtid)
    pointMap(el,'data/grt-20130215032034-%s/rays-%s.data'%(grtid,loc),'Antartica_average',title='Antartica (average speed distribution)')
    loc="lat_-7.78500e+01__lon_1.66400e+02"
    grtid="431628" #
    print "Probability:",np.loadtxt("data/grt-20130215032034-%s/geographic.prob"%grtid)
    pointMap(el,'data/grt-20130215032034-%s/rays-%s.data'%(grtid,loc),'Antartica_apex',title='Antartica (average speed distribution)')
    #"""

def genFireballsFile():
    #READ FIREBALL DATA
    spy.furnsh("kernels.txt")
    
    f=open('references/Data/fireballs.csv','r')
    lines=csv.reader(f,delimiter=',')
    i=0
    fields=\
    [
        'DateTime','TimeET',
        'Latitude','Longitude','Altitude', #Degrees
        'Velocity','vx','vy','vz', #km/s, with respect to J2000
        'TotalRadiatedEnergy', #Joules
        'TotalImpactEnergy' #ktons
    ]
    
    print "Reading csv data..."
    data=[]
    for line in lines:
        i+=1
        if i==1:continue
        j=-1
        subline=[]
        for value in line:
            j+=1
            field=fields[j]
            if value=='':
                subline+=[123456789]
                continue
            if field=='DateTime':
                et=spy.str2et(value)
                datetime=value.replace(" ","").replace(":","").replace("-","")
                date,time=value.split(" ")
                h,m,s=time.split(":")
                datetime=date.replace("-","")+"%02d%02d%02d"%(int(h),int(m),int(s))
                subline+=[int(datetime),et]
                j+=1
            elif field=="Latitude":
                sgn=+1
                if 'S' in value:sgn=-1
                value=sgn*float(value[:-1])
                subline+=[value]
            elif field=="Longitude":
                sgn=+1
                if 'W' in value:sgn=-1
                value=sgn*float(value[:-1])
                subline+=[value]
            else:
                subline+=[float(value)]
        data+=[subline]
    data=np.array(data)
    f.close()

    filename="util/data/fireballs.data"
    np.savetxt(filename,data,fmt="%-+26.17e")

    print "Creating header..."
    ffrm="#"
    fhead=()
    i=0
    for field in fields:
        ffrm+="%-26s "
        fhead+=("%d:"%i+field,)
        i+=1
    ffrm+="\n"
    f=open(filename,"a")
    f.write(ffrm%fhead)
    f.close()

def distributionFireballs(qload=0):
    data=np.loadtxt("util/data/fireballs.data")
    ovimps=data[:,5]
    lats=data[:,2]
    lons=data[:,3]
    energy=data[:,10]
    datetime=data[:,0]
    ets=data[:,1]

    condv=ovimps!=123456789
    vimps=ovimps[condv]

    print "Number of velocities:",len(vimps)

    conlatlon=(lats!=123456789)*(lons!=123456789)*(energy!=123456789)
    datetime=datetime[conlatlon]
    ets=ets[conlatlon]
    lats=lats[conlatlon]
    lons=lons[conlatlon]
    energy=energy[conlatlon]
    vimpqs=ovimps[conlatlon]

    print "Number of impacts:",len(lats)

    print "Most recent impact: %d"%datetime[0],lats[0],lons[0],energy[0]

    emin=energy.min()
    emax=energy.max()
    logemin=np.log10(emin)
    logemax=np.log10(emax)
    smin=2
    smax=10
    me=(smax-smin)/(logemax-logemin)
    print "Range of energies:",emin,emax

    # DIRECTION OF IMPACT SITE WITH RESPECT TO APEX
    if not qload:
        print "Computing apex distributions..."
        projs=[]
        for i in xrange(len(lats)):
            #if i<390:continue
            cmd="./whereisapex.exe ET %.9e %.6e %.6e > /dev/null"%(ets[i],lats[i],lons[i])
            out=System(cmd)
            latapex=float(out.split("\n")[6])
            lonapex=float(out.split("\n")[7])
            projs+=[[ets[i],lats[i],lons[i],latapex,lonapex]]
            """
            print "Command:",cmd
            print "Lat.apex, Lon.apex:",latapex,lonapex
            exit(0)
            """

        projs=np.array(projs)
        projs=np.savetxt("util/data/fireball-apex.data",projs)

    # APEX LATITUDE DISTRIBUTION
    print "Loading apex distributions..."
    projs=np.loadtxt("util/data/fireball-apex.data")
    latapexs=projs[:,3]
    colatapexs=90.0-latapexs
    print "Range of apex directions:",len(colatapexs),colatapexs.min(),colatapexs.max()

    h,q=np.histogram(colatapexs,15,normed=1)
    qm=(q[:-1]+q[1:])/2

    fig=plt.figure()
    ax=fig.gca()
    #ax.plot(90-qm,h/np.cos(np.pi/2-qm*DEG))
    hc=h
    qc=qm

    hc=h/np.cos((90-qm)*DEG)
    qc=90-qm
    
    bins,n=histOutline(hc,qm)
    ax.plot(bins,n)
    ax.set_xlim((0,180))
    fig.savefig(FIGDIR+"fireball-latapex-distribution.png")

    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # APEX CUMULATIVE LATITUDE DISTRIBUTION
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    Fcum=cumDistrib(qm,h)
    Fcum[-1,1]=1.0

    fig=plt.figure()
    ax=fig.gca()
    ax.plot(Fcum[:,0],Fcum[:,1],label=r'$N(<\theta_{\rm apex})$, Towards Antapex')
    ax.plot(180-Fcum[:,0],1-Fcum[:,1],label=r'$N(<180^\circ-\theta_{\rm apex})$, Towards Apex')
    ax.set_xlim((0,180))
    ax.legend(loc='lower rigth',fontsize=10)
    fig.savefig(FIGDIR+"fireball-latapex-cum-distribution.png")

    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # APEX LONGITUDE DISTRIBUTION
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    print "Loading apex distributions..."
    projs=np.loadtxt("util/data/fireball-apex.data")
    lonapexs=projs[:,4]
    print "Range of apex longitudes:",len(lonapexs),lonapexs.min(),lonapexs.max()

    h,q=np.histogram(lonapexs,30,normed=1)
    qm=(q[:-1]+q[1:])/2

    fig=plt.figure()
    ax=fig.gca()
    bins,n=histOutline(h,q)
    ax.plot(bins,n)
    ax.set_xlim((-180,180))
    fig.savefig(FIGDIR+"fireball-lonapex-distribution.png")

    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # MAP OF FIREBALLS AS A FUNCTION OF APEX DIRECTION
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    #"""
    fig=plt.figure()
    ax=fig.gca()
    proj='cyl'
    proj='moll'
    proj='vandg'
    proj='sinu'
    map=drawMap(proj=proj,proj_opts=dict(lon_0=0,ax=ax),
                pars=[-45,0,45],mers=[0,45,90,135,180,225,270,315],
                pars_opts=dict(labels=[0,1,0,0],fontsize=8),
                mers_opts=dict(labels=[0,0,0,1],fontsize=8))
    ms=5
    plotMap(map,lonapexs,latapexs,lw=0,
            marker='o',color='r',ms=ms,mec='none')

    fig.tight_layout()
    fig.savefig(FIGDIR+"fireballs-apex-map.png")
    #"""
    
    # IMPACT VELOCITY AS A FUNCTION OF TETAPEX
    fig=plt.figure()
    ax=fig.gca()
    ax.plot(colatapexs,vimpqs,'ko')

    vmax=vimpqs[vimpqs<123456789].max()
    vmin=vimpqs.min()
    ax.set_ylim((VESC,40))

    ax.set_xlabel(r"${\theta_{\rm apex}^{\rm geo}}^\circ$",fontsize=14)
    ax.set_ylabel(r"$v_{\rm imp}$ (km/s)",fontsize=14)
    fig.savefig(FIGDIR+"fireball-apex-velocity.png")

    exit(0)

    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    # MAP OF IMPACTS
    # &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    fig=plt.figure()
    ax=fig.gca()
    proj='moll'
    proj='robin'
    proj='cyl'
    map=drawMap(proj=proj,proj_opts=dict(lon_0=0,ax=ax),
                pars=[-45,0,45],mers=[0,45,90,135,180,225,270,315],
                pars_opts=dict(labels=[0,1,0,0],fontsize=8),
                mers_opts=dict(labels=[0,0,0,1],fontsize=8),
                coasts=True)
    for i in xrange(len(lats)):
        lenergy=np.log10(energy[i])
        ms=me*(lenergy-logemin)+smin
        plotMap(map,lons[i],lats[i],lw=0,
                marker='o',color='r',ms=ms,mec='none')

    fig.tight_layout()
    fig.savefig(FIGDIR+"fireballs-map.png")

    # DISTRIBUTION OF VELOCITIES
    fig=plt.figure()
    ax=fig.gca()
    ax.hist(vimps,10)
    fig.savefig(FIGDIR+"fireballs-velocities.png")

def apexAngleChelyabinskTunguska():

    #Tunguska
    out=System("./whereisapex.exe '06/30/1908 00:14:00 UTC' 60.917 101.95")
    cecl=float(out.split("\n")[2])
    capx=float(out.split("\n")[3])
    qapex=np.arccos(capx)*RAD
    print "Direction with respect to apex Tunguska event:",qapex

    #Chelyabinsk
    out=System("./whereisapex.exe '02/15/2013 03:02:34' 54.4 63.5")
    cecl=float(out.split("\n")[2])
    capx=float(out.split("\n")[3])
    qapex=np.arccos(capx)*RAD
    print "Direction with respect to apex Chelyabinsk event:",qapex

def testParticle():
    system("make && python throwaray.py 54.4 63.5 8e+04 17.664618 104.975030 -2.045864e+01 '02/15/2013 3:20:34 UTC' -0.5")
    system("cp scratch/ray-elements.png scratch/ray-orbit.png "+FIGDIR)

def mapProbability():

    #CHELYABINSK MAP
    #"""
    date="20130215032034"
    print "Mapping %s..."%date 
    #grtid="2B353A" #APEX
    grtid="3CAA5C" #AVERAGE
    qlat=54.4
    qlon=63.5
    if path.isfile("data/grt-%s-%s/Pmatrix.data"%(date,grtid)):qmatrix=0
    else:qmatrix=1
    cmd="python mapatsource.py data/grt-%s-%s %d %f %f"%(date,grtid,qmatrix,qlat,qlon)
    print "Excuting:",cmd
    system(cmd)
    system("cp data/grt-%s-%s/probability-map-contour.png %s/probability-map-contour-%s.png"%(date,grtid,FIGDIR,date))
    #"""

    #PRE CHELYABINSK MAP
    #"""
    date="20130214212034"
    print "Mapping %s..."%date 
    #grtid="9D35F9" #APEX
    grtid="C7BF6A" #AVERAGE
    qlat=54.4
    qlon=63.5
    if path.isfile("data/grt-%s-%s/Pmatrix.data"%(date,grtid)):qmatrix=0
    else:qmatrix=1
    cmd="python mapatsource.py data/grt-%s-%s %d %f %f"%(date,grtid,qmatrix,qlat,qlon)
    print "Excuting:",cmd
    system(cmd)
    system("cp data/grt-%s-%s/probability-map-contour.png %s/probability-map-contour-%s.png"%(date,grtid,FIGDIR,date))
    #"""

    #TUNGUSKA MAP
    #"""
    date="19080630001400"
    print "Mapping %s..."%date 
    #grtid="732A04" #APEX
    grtid="1BB07B" #AVERAGE
    qlat=60.917
    qlon=101.95
    if path.isfile("data/grt-%s-%s/Pmatrix.data"%(date,grtid)):qmatrix=0
    else:qmatrix=1
    cmd="python mapatsource.py data/grt-%s-%s %d %f %f"%(date,grtid,qmatrix,qlat,qlon)
    print "Excuting:",cmd
    system(cmd)
    system("cp data/grt-%s-%s/probability-map-contour.png %s/probability-map-contour-%s.png"%(date,grtid,FIGDIR,date))
    #"""
    
    #1963 EVENT
    #"""
    date="19630803164500"
    print "Mapping %s..."%date 
    #grtid="C535BA" #APEX
    grtid="A14D18" #AVERAGE
    qlat=-51.0
    qlon=+24.0
    if path.isfile("data/grt-%s-%s/Pmatrix.data"%(date,grtid)):qmatrix=0
    else:qmatrix=1
    cmd="python mapatsource.py data/grt-%s-%s %d %f %f"%(date,grtid,qmatrix,qlat,qlon)
    print "Excuting:",cmd
    system(cmd)
    system("cp data/grt-%s-%s/probability-map-contour.png %s/probability-map-contour-%s.png"%(date,grtid,FIGDIR,date))
    #"""

def probabilityFireballs():

    rundir="scratch/runs/"
    data=np.loadtxt("util/data/fireballs.data")
    lats=data[:,2]
    lons=data[:,3]
    datetime=data[:,0]
    ets=data[:,1]
    Es=data[:,10]
    vimps=data[:,5]

    cond=(lats!=123456789)*(lons!=123456789)

    lats=lats[cond]
    lons=lons[cond]
    datetime=datetime[cond]
    ets=ets[cond]
    Es=Es[cond]
    vimps=vimps[cond]
    
    nfire=len(lats)

    print "%d fireballs to analyse"%nfire

    nruns=4
    ngroup=nfire/nruns

    fs=[]
    for j in xrange(nruns):
        k=j+1
        fs+=[open(rundir+"run-%d.sh"%k,"w")]

    run=0
    j=0
    flog=open(rundir+"runs.log","w")
    for i in xrange(nfire):
        j+=1
        if (i%ngroup)==0 and run<nruns:
            run+=1
            j=1

        lat=lats[i]
        lon=lons[i]
        E=Es[i]
        if E==123456789:E=-1
        vimp=vimps[i]
        if vimp==123456789:vimp=-1

        date="%.0f"%datetime[i]
        Y=date[:4];x=4
        M=date[x:x+2];x+=2
        D=date[x:x+2];x+=2
        h=date[x:x+2];x+=2
        m=date[x:x+2];x+=2
        s=date[x:x+2];x+=2
        datestr="%s/%s/%s %s:%s:%s UTC"%(M,D,Y,h,m,s)

        QVEL=1
        NAME="Lat. %.4f, Lon. %.4f, vimp = %.2f, E = %.2f, Date %s"%(lat,lon,vimp,E,datestr)
        makestr="QVEL=%d & NAME=%s"%(QVEL,NAME)
        md5str=MD5STR(makestr,len=6)
        odir="data/grt-%s-%s"%(date,md5str)
        system("mkdir -p %s"%odir)

        fgeo=open(rundir+"geo.%d"%i,"w")
        fgeo.write("-1 -1 %.1f %.1f\n"%(lon,lat))
        fgeo.close()

        cmd="""echo;echo;echo '%s';echo 'Run %d/%d...';echo '%s';echo;echo
python makeagravray.py '%s' deg %s/geo.%d locals.dat %d '%s' > %s/grt.log
echo '%s' >> %s/runs.completed
"""%("*"*50,j,ngroup,"*"*50,datestr,rundir,i,QVEL,NAME,odir,date,rundir)
        
        fs[run-1].write(cmd+"\n")
        
        flog.write("%d %d %s\n"%(run,i,odir))
        #break

    for f in fs:f.close()
    flog.close()

def runProbabilityFireballs():

    rundir="scratch/runs/"
    data=np.loadtxt("util/data/fireballs.data")
    lats=data[:,2]
    lons=data[:,3]
    datetime=data[:,0]
    ets=data[:,1]
    Es=data[:,10]
    vimps=data[:,5]

    cond=(lats!=123456789)*(lons!=123456789)

    lats=lats[cond]
    lons=lons[cond]
    datetime=datetime[cond]
    ets=ets[cond]
    Es=Es[cond]
    vimps=vimps[cond]

    f=open("scratch/runs/runs.log")
    fp=open("scratch/fireball-probabilities.dat","w")
    for line in f:
        run,i,odir=line.strip("\n").split()
        i=int(i)
        """
        cmd="rm -r %s"%odir
        print cmd
        system(cmd)
        """
        # Read probability
        probs=np.loadtxt("%s/geographic.prob"%odir)
        fp.write("%d %.2f %.2f %.4f\n"%(i,probs[0],probs[1],probs[2]))
    fp.close()

    probs=np.loadtxt("scratch/fireball-probabilities.dat")

    fig=plt.figure()
    ax=fig.gca()
    ax.hist(probs[:,3],30)
    fig.savefig(FIGDIR+"fireball-probability-distribution.png")

def testProbability():
    
    from scipy import stats

    x=np.random.normal(loc=0,scale=1,size=1000)
    y=stats.norm.pdf(x)
    
    fig=plt.figure()
    ax=fig.gca()
    ax.hist(x,30,normed=1)
    ax.plot(x,y,'o',ms=0.5)
    fig.savefig(FIGDIR+"probability-test.png")

    fig=plt.figure()
    ax=fig.gca()
    ax.hist(y,30)
    fig.savefig(FIGDIR+"probability-pdf-test.png")

def plotConfigurationSpace(el):
    #LOAD NEOS DATA
    qes=el[:,0];qmin=qes.min();qmax=qes.max()
    ees=el[:,1];emin=ees.min();emax=ees.max()
    ies=np.log10(el[:,2]);imin=ies.min();imax=ies.max()
    ies=el[:,3];imin=ies.min();imax=ies.max()
    aes=el[:,4];amin=aes.min();amax=aes.max()
    qlow=qmin;qup=qmax
    elow=emin;eup=emax
    ilow=imin*0;iup=imax*0+0.5

    #CONFIGURATION SPACE
    cmap='rainbow'
    interpolation='nearest'
    bins=30
    factor=1.0
    pprop=dict(ms=8,mec='none')
    iprop=dict(ms=15,mec='k',color='w')
    fig=plt.figure(figsize=(18,7))
    fsize=18
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # q vs. e
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    condi=np.abs(qes)>=0
    H,xe,ye=np.histogram2d(qes[condi],ees[condi],bins=bins,normed=True)
    axae=fig.add_subplot(131)
    scale=(qup-qlow)/(eup-elow)
    img=axae.imshow(H.transpose(),origin='lower',
                    interpolation=interpolation,
                    extent=(qmin,qmax,emin,emax),aspect=scale/factor,cmap=cmap)
    axae.set_xlabel("$q$ (AU)",fontsize=fsize)
    axae.set_ylabel("$e$",fontsize=fsize)
    axae.set_xlim((qlow,qup))
    axae.set_ylim((elow,eup))
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # q vs. i
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    condi=np.abs(qes)>=0
    H,xe,ye=np.histogram2d(qes[condi],ies[condi],bins=bins,normed=True)
    axai=fig.add_subplot(132)
    scale=(qup-qlow)/(iup-ilow)
    img=axai.imshow(H.transpose(),origin='lower',
                    interpolation=interpolation,
                    extent=(qmin,qmax,imin,imax),aspect=scale/factor,cmap=cmap)
    axai.set_xlabel("$q$ (AU)",fontsize=fsize)
    axai.set_ylabel("$\log(i^\circ)$",fontsize=fsize)
    axai.set_xlim((qlow,qup))
    axai.set_ylim((ilow,iup))
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # e vs. i
    # $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    condi=np.abs(qes)>=0
    H,xe,ye=np.histogram2d(ees[condi],ies[condi],bins=bins,normed=True)
    axei=fig.add_subplot(133)
    scale=(eup-elow)/(iup-ilow)
    img=axei.imshow(H.transpose(),origin='lower',
                    interpolation=interpolation,
                    extent=(emin,emax,imin,imax),aspect=scale/factor,cmap=cmap)
    axei.set_xlabel("$e$",fontsize=fsize)
    axei.set_ylabel("$\sin(i^\circ)$",fontsize=fsize)
    axei.set_xlim((elow,eup))
    axei.set_ylim((ilow,iup))
    return axae,axai,axei

def visualizeProcess(el):

    """
    Location of the antapex: lat. 19.5, lon. -141.4
    """

    #INPUT INFO
    date="02/15/2013 03:20:34 UTC"
    out=System("./whattimeisit.exe '%s' ET > /dev/null"%date)
    t=float(out.split("\n")[4])
    alt=8e4

    #HAWAII
    odir="data/grt-20130215032034-ED2863/";lat=+21.8;lon=-157.0
    #CHELYABINSK
    #odir="data/grt-20130215032034-FA0C86/";lat=+54.4;lon=+63.5
    #MADAGASCAR
    #odir="data/grt-20130215032034-21A94B/";lat=-18.9;lon=+47.5

    #PARAMETERS
    #dmax=0.1
    dmax=0.15
    sigma=wNormalization(dmax)
    wmax=sigma*wFunction(0,dmax)
    normal=2000.0
    fparam=(0.9721768,6.84870896,2.40674371)

    #READ DATA
    sites=np.loadtxt(odir+"locals.dat")

    k=0
    hold=0
    caz=0
    qps=[];eps=[];ips=[]
    for site in sites:
        plt.close("all")

        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        #Info about direction
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        h,A,v0,v1,v2,v3,end=site

        if h==hold:
            if caz==3:continue
        else:
            if caz==3:caz=0
        hold=h
        caz+=1

        print "*"*80,"\nSite %d (caz = %d)\n"%(k,caz),"*"*80

        print "Direction:",site
        print "Zenith angle:",90-h
        f=open("tmp/locals.data","w")
        f.write("#\n");

        #Test
        #h=18.0;A=105.0;v0=v1=v2=v3=20.0;
        f.write("%-+20.4e%-+20.4e%-+20.4e%-+20.4e%-+20.4e%-+20.4e%-+20d\n"%(h,A,v0,v1,v2,v3,end))
        f.close()

        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        #1-Direction in the sky
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        """
        proj='moll'
        fig=plt.figure();ax=fig.gca()
        map_dir=drawMap(proj=proj,proj_opts=dict(lon_0=180,ax=ax),
                        pars=[-45,0,45],mers=[0,45,90,135,180,225,270,315],
                        pars_opts=dict(labels=[1,1,0,0],fontsize=8),
                        mers_opts=dict(labels=[0,0,0,0],fontsize=8))
        plotMap(map_dir,A,h,lw=0,
                marker='o',color='r',ms=4,mec='none')
        plt.show()
        """
        
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        #2-Apex direction
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        cmd="./whereisapex.exe '%s' %f %f 90.0 0.0 20.0 > /dev/null"%(date,lat,lon)
        print "Executing:",cmd
        out=out2dict(System(cmd))
        print "Geographic Lat. apex = ",90-float(out["LATAPEX"])
        print "Velocity (20 km/s) Lat. apex = ",90-float(out["LATAPEX"])

        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        #3-Integration
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        QVEL=1
        inifile="tmp/locals.data"
        outfile="tmp/rays.data"
        cmd="./throwrays.exe %.9e %.5e %.5e %.4e %s %d %s"%(t,lat,lon,alt,inifile,QVEL,outfile)
        print "Executing:",cmd
        system(cmd)
        ray=np.loadtxt(outfile)
        try:
            h,Az,vimp,xobs,yobs,zobs,vxobs,vyobs,vzobs,q,e,i,Omega,omega,M,qapex=ray
        except:
            continue

        a=q/(1-e)
        print "Orbital elements:",q,a,e,i,Omega,omega
        print "Direction with respect to apex:",qapex
        flux=theoFlux_DoubleTrigCos(qapex,*fparam)
        print "Flux of objects in that direction:",flux
        
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        #4-Density
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        verb=0
        #distform=drummondDistance(q,e,i)
        distform=zappalaDistance(a,e,np.sin(i*DEG),Omega,omega)
        result=np.array(mysqlSelect("%s, Perihelion_dist, e, i, sini, a, Node, Peri"%distform,
                                    "NEOS",
                                    "where %s<%e order by %s"%(distform,(2*dmax)**2,distform),"array"))
        ntarg=result.shape[0]
        print "Number of close objects:",ntarg
        """
        for target in result:
            d2,qt,et,it,sinit,at,Ot,ot=target
            d=d2**0.5
            print "q=%.3f,%.3f"%(q,qt),"e=%.3f,%.3f"%(e,et),"i=%.3f,%.3f"%(i,it),"sini=%.3f,%.3f"%(np.sin(i*DEG),sinit),"a=%.3f,%.3f"%(a,at),"O = %.3f,%.3f"%(Omega,Ot),"o = %.3f,%.3f"%(omega,ot),"d = %.3f"%d
            ka=5./4
            ke=ki=2
            kw=kO=1e-4
            am=(a+at)/2
            d2c=1/np.sqrt(am)*(ka*((at-a)/am)**2+ke*(et-e)**2+ki*(sinit-np.sin(i*DEG))**2+kO*(Omega-Ot)**2+kw*(omega-ot)**2)
            dc=np.sqrt(d2c)
            print d2c,dc
        return
        """

        density=0
        if ntarg>0:
            n=0
            for target in result:
                d2,qt,et,it,sinit,at,Ot,ot=target
                d=d2**0.5
                p=sigma*wFunction(d,dmax)
                print "q=%.3f,%.3f"%(q,qt),"e=%.3f,%.3f"%(e,et),"i=%.3f,%.3f"%(i,it),"sini=%.3f,%.3f"%(np.sin(i*DEG),sinit),"a=%.3f,%.3f"%(a,at),"O = %.3f,%.3f"%(Omega,Ot),"o = %.3f,%.3f"%(omega,ot),"d = %.3f"%d,"p=",p
                density+=p
                n+=1
            print "Density:",density
        else:
            density=0

        #7-Probability
        Pu=density/wmax
        Pn=flux*Pu

        print "Probability without flux correction:",Pu
        print "Probability with flux correction:",Pn

        """
        cmd="python throwaray.py %f %f %f %f %f -%f '%s' -0.50 100 1"%(lat,lon,alt,h,Az,vimp,date)
        print "Executing:",cmd
        system(cmd)
        """

        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        #5-Position in configuration & physical space
        #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        #GENERATE DENSITY MAP
        funi=lambda i:np.sin(i*np.pi/180)
        axae,axai,axei=plotConfigurationSpace(el)

        axae.plot(q,e,'ko',ms=10)
        axai.plot(q,funi(i),'ko',ms=10)
        axei.plot(e,funi(i),'ko',ms=10)

        axae.plot(qps,eps,'kv',ms=5)
        axai.plot(qps,ips,'kv',ms=5)
        axei.plot(eps,ips,'kv',ms=5)

        qps+=[q];eps+=[e];ips+=[np.log10(i)]

        axae.text(0.9*q,0.9*e,"%.1f,%.1f,%.1f,%.1f"%(qapex,flux,Pu,Pn),fontsize=8,ha='right')

        plt.show()

        #raw_input()
        #if k>5:break
        #break
        k+=1

def testUniformOmega(el):
    qes=el[:,0];
    ees=el[:,1];
    ies=np.log10(el[:,2]);
    aes=qes/(1-ees**2);
    Omegas=el[:,3]
    omegas=el[:,4]

    #"""
    #cond=(qes<=0.7)*(qes>0.6)
    cond=(qes>0.8)*(ees<0.6)*(ees>0.4)*(ies>0.5)*(ies<1)
    Omegas=Omegas[cond]
    omegas=omegas[cond]
    #"""
    nobs=len(Omegas)

    hO,O=np.histogram(Omegas,50,normed=1)
    ho,o=np.histogram(omegas,50,normed=1)
    
    fig=plt.figure()
    ax=fig.gca()
    bins,n=histOutline(hO,O)
    ax.plot(bins,n,label=r'$\Omega$ Distribution')
    bins,n=histOutline(ho,o)
    ax.plot(bins,n/2,label=r'$\omega$ Distribution')
    ax.set_xlim((0,360))
    ax.set_title('%d EC-NEOs (Sept.2016)'%nobs,position=(0.5,1.02))
    ax.set_xlabel('Angle ($^\circ$)')
    ax.set_ylabel('Frequency')
    ax.set_yticks([])
    ax.legend(loc='lower right')
    fig.tight_layout()
    fig.savefig(FIGDIR+"OrientationDistribution.png")

#############################################################
#EXECUTE
#############################################################
#distanceTunguskaChelyabisnk()
#plotGeographicPositions()
#experimen2(elements)
