from gravray import *

qload=1
try:qload=int(argv[1])
except:pass

if qload>0:
    props="Perihelion_dist,Aphelion_dist,e,i"
    condition="where (Perihelion_dist<=0.9833 and Aphelion_dist>=1.0167) or (Perihelion_dist<=1.0167 and Perihelion_dist>=0.9833) or (Aphelion_dist<=1.0167 and Aphelion_dist>=0.9833) or (Aphelion_dist<=1.0167 and Perihelion_dist>=0.9833) and (e>0 and e<1)"
    listdict=mysqlSelect(props,
                         condition=condition)
    elements=listdict2matrix(listdict,keys=props.split(","))

def velocityDistribution(el):

    verbose=0

    nbodies=el.shape[0]

    #OBRITAL VELOCITY
    vorb=np.sqrt(6.67E-11*2E30/1.5E11)/1E3
    print "Earth orbital velocity (km/s): ",vorb

    #ESCAPE VELOCITY
    vesc2=(2*5.98E24/2E30)*(1.5e8/6.34e3)
    vesc=np.sqrt(vesc2)*vorb
    print "Earth escape velocity (v^2,v,v[km/s]): ",vesc2,np.sqrt(vesc2),vesc

    allindexes=np.arange(nbodies)
    nsamples=30
    nsubsample=3000
    nbins=30
    hmean=np.zeros(nbins)

    for n in xrange(nsamples):
        indexes=np.random.choice(allindexes,nsubsample)
        velocities=[]
        for i in indexes:

            #ORBITAL ELEMENTS
            q=el[i,0]
            Q=el[i,1]
            e=el[i,2]
            i=el[i,3]

            a=q/(1-e)
            if verbose:print "Orbital elements (a,e,i):",a,e,i

            #ASTEROID VELOCITY
            r=1.0
            va2=(2/r-1/a)
            if verbose:print "va^2 = ",va2

            #HORIZONTAL COMPONENT VELOCITY
            h=np.sqrt(a*(1-e**2))
            ydot=h*np.cos(i*DEG)/r
            if verbose:print "ydot = ",ydot

            #V-INFINITY
            vinf2=va2+1-2*ydot
            if verbose:print "vinf^2 = ",vinf2

            #V-IMP
            vimp2=vinf2+vesc2
            if verbose:print "vimp^2 = ",vimp2
            vimp=np.sqrt(vimp2)
            if verbose:print "vimp = ",vimp

            if verbose:print "vimp [km/s] = ",vimp*vorb
            velocities+=[vimp*vorb]

        velocities=np.array(velocities)
        h,x=np.histogram(velocities,nbins,(vesc,4*vesc))
        hmean+=h

    ###################################################
    #SAVE CUMMULATIVE IMPACT VELOCITY DISTRIBUTION
    ###################################################
    xbin=x[:-1]
    hmean/=1.0*nsamples
    Fmean=np.zeros_like(hmean)
    for i in xrange(1,hmean.shape[0]):
        Fmean[i]=Fmean[i-1]+hmean[i-1]
    Fmean/=nsubsample
    np.savetxt("velocity-cum.data",np.vstack(
        (np.append(xbin,x[-1]),
         np.append(Fmean,[1]))
    ).transpose())

    ###################################################
    #SAVE HISTOGRAM
    ###################################################
    np.savetxt("velocity-pdf.data",np.vstack((xbin,hmean/(1.0*nsubsample))).transpose())

    ###################################################
    #PLOTS
    ###################################################
    fig=plt.figure()
    ax=fig.gca()
    ax.plot(xbin,hmean)
    ax.axvline(vesc);
    fig.savefig("scratch/velocities-pdf.png")

    fig=plt.figure()
    ax=fig.gca()
    ax.plot(xbin,Fmean)
    ax.axvline(vesc);
    fig.savefig("scratch/velocities-cum.png")

def generateVelocitiesPlot():
    nsample=10000
    v=generateVelocities("velocity-cum.data",nsample)

    velpdf=np.loadtxt("velocity-pdf.data")
    nbins=velpdf.shape[0]
    vmin=velpdf[:,0].min();vmax=velpdf[:,0].max();

    h,x=np.histogram(v,nbins,(vmin,vmax))
    fig=plt.figure()
    ax=fig.gca()
    ax.plot(x[:-1],h,'+')
    ax.plot(velpdf[:,0],velpdf[:,1]*nsample)
    fig.savefig("scratch/velocities-pdf-theo.png")
