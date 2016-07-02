from gravray import *

#############################################################
#INPUTS
#############################################################
#Ensamble directory
edir=argv[1]

#Load data?
qload=1
try:qload=int(argv[2])
except:pass

#Load matrix?
qmat=1
try:qmat=int(argv[3])
except:pass

if qload:
    #############################################################
    #GET ENSAMBLE PROPERTIES
    #############################################################
    out=System("ls -m %s/*.data.prob"%edir)
    dataensamble=[]
    for fname in out.split(","):
        fname=fname.strip("\n ");
        pcoords=fname.split("-lat")[1].split("__")
        lat=float(pcoords[0].replace("_",""))
        lon=float(pcoords[1].replace("lon_","").replace(".data.prob",""))
        data=np.loadtxt(fname)
        Ptot=data[:,7].sum()
        dataensamble+=[[lat,lon,Ptot]]
    dataensamble=np.array(dataensamble)

def readEnsamble(de):

    #READ OBSERVER MATRICES
    obsdata=np.loadtxt("%s/observers-matrices.dat"%edir)
    latsize=int(obsdata[:,0].max())+1
    lonsize=int(obsdata[:,1].max())+1
    mlats=np.unique(obsdata[:,2])
    mlons=np.unique(obsdata[:,3])
    colat=np.zeros((latsize,lonsize))
    coapex=np.zeros((latsize,lonsize))
    k=0
    for i in xrange(latsize):
        for j in xrange(lonsize):
            colat[i,j]=obsdata[k,4]
            coapex[i,j]=obsdata[k,5]
            #print obsdata[k,2],obsdata[k,3],colat[i,j],coapex[i,j]
            #raw_input()
            k+=1

    #############################################################
    #GENERAL PROPERTIES
    #############################################################
    lats=de[:,0]
    lons=de[:,1]
    Ptots=de[:,2]
    npoints=len(Ptots)
    
    Pmin=Ptots.min()
    Pmax=Ptots.max()
    print "Ilumination range:",Pmin,Pmax

    #############################################################
    #CONTOUR
    #############################################################
    dlat=5.0
    dlon=5.0
    glats=np.arange(-90.,90.+dlat,dlat)
    glons=np.arange(0.,360.+dlat,dlat)

    Nlats=len(glats)
    Nlons=len(glons)
    Pmatrix=np.zeros((Nlats,Nlons))
    lonms=np.mod(lons,360.0)

    if qmat:
        fm=open("%s/Pmatrix-human.data"%edir,"w")
        for i in xrange(Nlats):
            for j in xrange(Nlons):
                glat=glats[i]
                glon=glons[j]
                alphas=np.array([arcDistance([glon*DEG,glat*DEG],[lon*DEG,lat*DEG])*RAD for lon,lat in zip(lonms,lats)])
                n=alphas.argsort()[0]
                Pmatrix[i,j]=Ptots[n]
                fm.write("%-06d%-06d%-+12.5f%-+12.5f%-+12.3e\n"%(i,j,glat,glon,Ptots[n]))
        fm.close()
        np.savetxt("%s/Pmatrix.data"%edir,Pmatrix)
    else:
        Pmatrix=np.loadtxt("%s/Pmatrix.data"%edir)

    #############################################################
    #PLOT DEPENDENCE ON LONGITUDE
    #############################################################
    fig=plt.figure()
    ax=fig.gca()
    ax.plot(lons,Ptots,'ro')
    ax.set_xlabel("Longitude (degrees)")
    ax.set_ylabel("Normalized illumination")
    fig.savefig("%s/ilumination-longitudes.png"%edir)

    #############################################################
    #PLOT MAP
    #############################################################
    plt.close("all")
    proj='robin'
    map=drawMap(proj=proj,proj_opts=dict(lon_0=0),
                pars=[-45,0,45],mers=[0,45,90,135,180,225,270,315],
                pars_opts=dict(labels=[1,1,0,0],fontsize=8),
                mers_opts=dict(labels=[0,0,0,1],fontsize=8),
                coasts=True
            )
    ax=plt.gca()

    colormap="spectral"
    cmap=cm.get_cmap(colormap)
    for i in xrange(npoints):
        lon=np.mod(lons[i],360.0)
        lat=lats[i]
        P=Ptots[i]
        Pnorm=(P-Pmin)/(Pmax-Pmin)
        color=cmap(Pnorm)
        x,y=map(lon,lat)
        plt.plot(x,y,'o',color=color,mec='k')
    plt.savefig("%s/illumination-map-discrete.png"%edir)


    #############################################################
    #PLOT CONTOUR OF ECLIPTIC COLATITUDE
    #############################################################
    plt.close("all")
    fig=plt.figure()
    proj='robin'
    map=drawMap(proj=proj,proj_opts=dict(lon_0=0.0),
                pars=[-45,0,45],mers=[0,45,90,135,180,225,270,315],
                pars_opts=dict(labels=[1,1,0,0],fontsize=8),
                mers_opts=dict(labels=[0,0,0,1],fontsize=8),
                coasts=True
            )
    ax=plt.gca()

    #SHIFT P-MATRIX
    Pmatrix,glonds=Grid(180.,Pmatrix,glons,start=False)
    LONS,LATS=np.meshgrid(glonds,glats)
    LN,LT=map(LONS,LATS)

    #SHIFT COLATITUDE
    colat,mlonds=Grid(180.,colat,mlons,start=False)
    LONS,LATS=np.meshgrid(mlonds,mlats)
    cLN,cLT=map(LONS,LATS)

    #SHIFT COAPEX
    coapex,mlonds=Grid(180.,coapex,mlons,start=False)
    LONS,LATS=np.meshgrid(mlonds,mlats)
    aLN,aLT=map(LONS,LATS)

    #==============================
    #CONTOURS OF P
    #==============================
    #colormap="spectral"
    #colormap="gray"
    #colormap="jet"
    #colormap="rainbow"
    colormap="seismic"
    cmap=cm.get_cmap(colormap)
    levels=np.linspace(Pmin,Pmax,10)
    c=map.contourf(LN,LT,Pmatrix,levels=levels,alpha=0.5,lw=0,cmap=cmap)

    #==============================
    #CONTOUR OF COLATITUDES
    #==============================
    map.contour(cLN,cLT,colat,levels=[0],colors=['k'],linewidths=[3])

    #==============================
    #CONTOUR OF COAPEX
    #==============================
    map.contour(aLN,aLT,coapex,levels=[0,0.999,-0.999],colors=['k','r','b'],linestyles=['--','-','-'],linewidths=[2])

    #==============================
    #DECORATION
    #==============================
    plt.title("%s"%edir,position=(0.5,1.05))
    cax=plt.axes([0.05,0.1,0.9,0.05])
    cbar=plt.colorbar(c,drawedges=False,cax=cax,orientation='horizontal',
                      format='%.2f')
    cbar.ax.tick_params(labelsize=8)
    cbar.ax.set_title("Normalized illumination",fontsize=8)
    plt.savefig("%s/illumination-map-contour.png"%edir)

if qload:
    readEnsamble(dataensamble)

