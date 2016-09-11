from gravray import *

qload=1
try:qload=int(argv[1])
except:pass

if qload>0:
    props="Perihelion_dist, Aphelion_dist, e, i"

    """
    condition="where (Perihelion_dist<=0.9833 and Aphelion_dist>=1.0167) or (Perihelion_dist<=1.0167 and Perihelion_dist>=0.9833) or (Aphelion_dist<=1.0167 and Aphelion_dist>=0.9833) or (Aphelion_dist<=1.0167 and Perihelion_dist>=0.9833) and (e>0 and e<1)"
    """
    condition="where Perihelion_dist<=1 and Perihelion_dist>=(1-e)/(1+e) and (e>0 and e<1)"

    listdict=mysqlSelect(props,
                         condition=condition)
    elements=listdict2matrix(listdict,keys=props.split(", "))

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

    #SUN ESCAPE VELOCITY
    vescsun=np.sqrt((2*6.67e-11*2E30)/1.5e11)/1e3
    print "Sun escape velocity (v[km/s]): ",vescsun

    allindexes=np.arange(nbodies)
    nsamples=30
    nsubsample=6000
    nbins=60

    hmean=np.zeros(nbins)
    hmeaninfs=np.zeros(nbins)
    hmeanasts=np.zeros(nbins)
    for n in xrange(nsamples):
        indexes=np.random.choice(allindexes,nsubsample)
        velocities=[]
        velinfs=[]
        velasts=[]
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
            va=np.sqrt(va2)
            velasts+=[va*vorb]

            #HORIZONTAL Y-COMPONENT OF VELOCITY
            h=np.sqrt(a*(1-e**2))
            ydot=h*np.cos(i*DEG)/r
            if verbose:print "ydot = ",ydot

            #V-INFINITY
            vinf2=va2+1-2*ydot
            if verbose:print "vinf^2 = ",vinf2
            if vinf2>0:vinf=np.sqrt(vinf2)
            else:vinf=0
            velinfs+=[vinf*vorb]

            #V-IMP
            vimp2=vinf2+vesc2
            if verbose:print "vimp^2 = ",vimp2
            vimp=np.sqrt(vimp2)
            if verbose:print "vimp = ",vimp

            if verbose:print "vimp [km/s] = ",vimp*vorb
            velocities+=[vimp*vorb]

        velocities=np.array(velocities)
        velinfs=np.array(velinfs)

        h,x=np.histogram(velocities,nbins,(vesc,4*vesc))
        hmean+=h

        h,xinf=np.histogram(velinfs,nbins,(0,4*vesc))
        hmeaninfs+=h

        h,xast=np.histogram(velasts,nbins,(0.5*vorb,vescsun))
        hmeanasts+=h

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
    #SAVE CUMMULATIVE VELOCITY AT INFINITE DISTRIBUTION
    ###################################################
    xbin=xinf[:-1]
    hmeaninfs/=1.0*nsamples
    Fmean=np.zeros_like(hmeaninfs)
    for i in xrange(1,hmeaninfs.shape[0]):
        Fmean[i]=Fmean[i-1]+hmeaninfs[i-1]
    Fmean/=nsubsample
    np.savetxt("velinf-cum.data",np.vstack(
        (np.append(xbin,xinf[-1]),
         np.append(Fmean,[1]))
    ).transpose())

    ###################################################
    #SAVE HISTOGRAM VELOCITY AT INFINITE
    ###################################################
    np.savetxt("velinf-pdf.data",np.vstack((xbin,hmeaninfs/(1.0*nsubsample))).transpose())

    ###################################################
    #SAVE CUMMULATIVE ASTEROID VELOCITY DISTRIBUTION
    ###################################################
    xbin=xast[:-1]
    hmeanasts/=1.0*nsamples
    Fmean=np.zeros_like(hmeanasts)
    for i in xrange(1,hmeanasts.shape[0]):
        Fmean[i]=Fmean[i-1]+hmeanasts[i-1]
    Fmean/=nsubsample
    np.savetxt("velast-cum.data",np.vstack(
        (np.append(xbin,xast[-1]),
         np.append(Fmean,[1]))
    ).transpose())

    ###################################################
    #SAVE HISTOGRAM ASTEROID VELOCITY
    ###################################################
    np.savetxt("velast-pdf.data",np.vstack((xbin,hmeanasts/(1.0*nsubsample))).transpose())

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

def velocityDistributionApex(el):

    verbose=1
    nbodies=el.shape[0]

    #OBRITAL VELOCITY
    vorb=np.sqrt(6.67E-11*2E30/1.5E11)/1E3
    print "Earth orbital velocity (km/s): ",vorb

    #ESCAPE VELOCITY
    vesc2=(2*5.98E24/2E30)*(1.5e8/6.34e3)
    vesc=np.sqrt(vesc2)*vorb
    print "Earth escape velocity (v^2,v,v[km/s]): ",vesc2,np.sqrt(vesc2),vesc

    #SUN ESCAPE VELOCITY
    vescsun=np.sqrt((2*6.67e-11*2E30)/1.5e11)/1e3
    print "Sun escape velocity (v[km/s]): ",vescsun

    allindexes=np.arange(nbodies)

    #nsamples=30
    #nsubsample=6000

    nsamples=1
    nsubsample=1

    ngrid=10
    cvals=np.linspace(-1,1,ngrid+1)
    dgrid=cvals[1]-cvals[0]

    #Statistics matrix
    nvels=100
    Mapex=np.zeros((ngrid,ngrid,nvels))
    Napex=np.zeros((ngrid,ngrid))

    #Generating samples
    for n in xrange(nsamples):

        indexes=np.random.choice(allindexes,nsubsample)

        for i in indexes:

            #ORBITAL ELEMENTS
            q=el[i,0]
            Q=el[i,1]
            e=el[i,2]
            i=el[i,3]*DEG

            #e=0;q=1;Q=q*(1+e)/(1-e);i=0
            #e=0.2;q=0.9;Q=q*(1+e)/(1-e);i=30*DEG

            p=q*(1+e)
            hop=1/np.sqrt(p);
            a=q/(1-e)

            #ASCENDING OR DESCENDING NODE
            O=0
            if np.random.rand()>0.5:O=180*DEG
            
            #O=0
            #O=180*DEG

            #PERIHELIA ARGUMENT
            if O==0:
                w=np.arccos((p-1)/e)
                wpf=0.0*DEG
            else:
                w=np.arccos((1-p)/e)
                wpf=180*DEG

            if verbose:
                print "Mock asteroid:"
                print "\t","q,Q,a,p,e,i = ",q,Q,a,p,e,i*RAD
                print "\t","Omega = ",O*RAD
                print "\t","omega = ",w*RAD

            #ASTEROID VELOCITY (VIS VIVA)
            r=1.0
            va2=(2/r-1/a)
            va=np.sqrt(va2)

            #VELOCITY COMPONENTS
            xdot=-hop*e*np.cos(O)*np.sin(w)
            ydot=+hop*np.cos(O)*np.cos(i)*(np.cos(wpf)+e*np.cos(w))
            zdot=+hop*np.sin(i)*(np.cos(wpf)+e*np.cos(w))

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
                print "\t","Relative speed (velocity at infinite) = ",vinf2

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
            yapex=-zdotrel
            zapex=+ydotrel

            tapex=np.arccos(zapex/vinf)
            fapex=np.arctan2(yapex,xapex)

            costapex=zapex/vinf
            nfapex=fapex/PI

            print cvals
            costapex=0.9

            iapex=int((costapex+1)/dgrid)
            japex=int((nfapex+1)/dgrid)

            if verbose:
                print "Apex direction:"
                print "\t","xapex,yapex,zapex = ",xapex,yapex,zapex
                print "\t","Theta apex = ",tapex*RAD
                print "\t","Phi tapex = ",fapex*RAD
                print "\t","cos(Theta apex) = ",costapex
                print "\t","Phi tapex / Pi = ",nfapex
                print "\t","i,j = ",iapex,japex

            exit(0)
            velocities+=[vimp*vorb]

        velocities=np.array(velocities)
        velinfs=np.array(velinfs)

        h,x=np.histogram(velocities,nbins,(vesc,4*vesc))
        hmean+=h

        h,xinf=np.histogram(velinfs,nbins,(0,4*vesc))
        hmeaninfs+=h

        h,xast=np.histogram(velasts,nbins,(0.5*vorb,vescsun))
        hmeanasts+=h

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
    #SAVE CUMMULATIVE VELOCITY AT INFINITE DISTRIBUTION
    ###################################################
    xbin=xinf[:-1]
    hmeaninfs/=1.0*nsamples
    Fmean=np.zeros_like(hmeaninfs)
    for i in xrange(1,hmeaninfs.shape[0]):
        Fmean[i]=Fmean[i-1]+hmeaninfs[i-1]
    Fmean/=nsubsample
    np.savetxt("velinf-cum.data",np.vstack(
        (np.append(xbin,xinf[-1]),
         np.append(Fmean,[1]))
    ).transpose())

    ###################################################
    #SAVE HISTOGRAM VELOCITY AT INFINITE
    ###################################################
    np.savetxt("velinf-pdf.data",np.vstack((xbin,hmeaninfs/(1.0*nsubsample))).transpose())

    ###################################################
    #SAVE CUMMULATIVE ASTEROID VELOCITY DISTRIBUTION
    ###################################################
    xbin=xast[:-1]
    hmeanasts/=1.0*nsamples
    Fmean=np.zeros_like(hmeanasts)
    for i in xrange(1,hmeanasts.shape[0]):
        Fmean[i]=Fmean[i-1]+hmeanasts[i-1]
    Fmean/=nsubsample
    np.savetxt("velast-cum.data",np.vstack(
        (np.append(xbin,xast[-1]),
         np.append(Fmean,[1]))
    ).transpose())

    ###################################################
    #SAVE HISTOGRAM ASTEROID VELOCITY
    ###################################################
    np.savetxt("velast-pdf.data",np.vstack((xbin,hmeanasts/(1.0*nsubsample))).transpose())

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

def velocityDistributionFit(elements):
    pass

#velocityDistribution(elements)
