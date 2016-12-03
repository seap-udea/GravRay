from gravray import *

#############################################################
#INPUTS
#############################################################
#Ensamble directory
iarg=1
edir=argv[iarg];iarg+=1
inidata=np.loadtxt("%s/locals.dat"%edir)
Ninitial=len(inidata)

#Calculate matrix?
qmat=1
try:qmat=int(argv[iarg]);iarg+=1
except:pass

qspecial=0
try:
    qspecial=1
    qlat=float(argv[iarg]);iarg+=1
    qlon=float(argv[iarg]);iarg+=1
except:pass    

#############################################################
#LOAD PROBABILITIES
#############################################################
geofile=System("ls %s/geographic*.prob"%edir)
de=np.loadtxt(geofile)

#############################################################
#READ OBSERVER MATRICES
#############################################################
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
print "Plotting impact probability in latitude and longitude"
fig=plt.figure()
ax=fig.gca()
ax.plot(lons,Ptots,'ro')
ax.set_xlabel("Longitude (degrees)")
ax.set_ylabel("Impact probability")
fig.savefig("%s/probability-longitudes.png"%edir)

fig=plt.figure()
ax=fig.gca()
ax.plot(lats,Ptots,'ro')
ax.set_xlabel("Latitude (degrees)")
ax.set_ylabel("Impact probability")
fig.savefig("%s/probability-latitudes.png"%edir)

#############################################################
#PLOT MAP
#############################################################
print "Plotting impact probability map (discrete)"
plt.close("all")

fig=plt.figure(figsize=(8,5))
proj='robin'
map=drawMap(proj=proj,proj_opts=dict(lon_0=0),
            pars=[-45,-30,-15,0,15,30,45],mers=[0,45,90,135,180,225,270,315],
            pars_opts=dict(labels=[0,0,0,0],fontsize=8),
            mers_opts=dict(labels=[0,0,0,0],fontsize=8),
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

if qspecial:
    print "Plotting special location: (%f,%f)"%(qlat,qlon)
    x,y=map(qlon,qlat)
    ax.plot(x,y,'k*',ms=10)

date=edir.split("/")[1].split("-")[1]
ax.set_title("Impact probabilities at %s"%date,position=(0.5,1.05))

fig.tight_layout()
plt.savefig("%s/probability-map-discrete.png"%edir)

#############################################################
#PLOT CONTOUR OF ECLIPTIC COLATITUDE
#############################################################
print "Plotting impact probability map (continuous)"

plt.close("all")
fig=plt.figure(figsize=(8,6))
proj='robin'
pars=np.arange(-75,90,15)
mers=np.arange(0,360,30)
map=drawMap(proj=proj,proj_opts=dict(lon_0=0.0),
            pars=pars,mers=mers,
            pars_opts=dict(labels=[0,0,1,0],fontsize=8),
            mers_opts=dict(labels=[0,0,0,0],fontsize=8),
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
colormap="rainbow"
#colormap="seismic"
#colormap="Paired"
#colormap="PuOr"
#colormap="gnuplot"
#colormap="RdYlGn"
cmap=cm.get_cmap(colormap)
levels=np.linspace(Pmin,Pmax,1000)
#levels=np.linspace(0.3,0.7,1000)
c=map.contourf(LN,LT,Pmatrix,levels=levels,alpha=1.0,lw=0,cmap=cmap)

#==============================
#CONTOUR OF COLATITUDES
#==============================
map.contour(cLN,cLT,colat,levels=[0],colors=['k'],linewidths=[3])

#==============================
#CONTOUR OF COAPEX
#==============================
map.contour(aLN,aLT,coapex,levels=[0,0.999,-0.999],colors=['k','r','b'],linestyles=['--','-','-'],linewidths=[2])
plt.tight_layout()

#==============================
#SPECIAL LOCATION POINT
#==============================
if qspecial:
    print "Plotting special location: (%f,%f)"%(qlat,qlon)
    x,y=map(qlon,qlat)
    ax.plot(x,y,'k*',ms=10)
    
    alphas=np.array([arcDistance([qlon*DEG,qlat*DEG],[lon*DEG,lat*DEG])*RAD for lon,lat in zip(lonms,lats)])
    n=alphas.argsort()[0]
    Pspecial=Ptots[n]
    print "Probability of special location:",Pspecial

#==============================
#DECORATION
#==============================
date=edir.split("/")[1].split("-")[1]
ax.set_title("Impact probabilities at %s"%date,position=(0.5,1.05))
cax=plt.axes([0.05,0.1,0.9,0.05])
cbar=plt.colorbar(c,drawedges=False,cax=cax,orientation='horizontal',
                  format='%.2f')
cbar.ax.tick_params(labelsize=8)
#cbar.ax.axvline(Pspecial,lw=3,color='k')
cbar.ax.plot((Pspecial-Pmin)/(Pmax-Pmin),0.5,'k*',ms=10)
cbar.ax.set_title("Normalized probability",fontsize=14,position=(0.5,-1.5))

#fig.tight_layout()
plt.savefig("%s/probability-map-contour.png"%edir)

#############################################################
#REMOVE DATA
#############################################################
del obsdata
del Pmatrix
