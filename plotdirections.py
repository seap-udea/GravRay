from gravray import *

###################################################
#PARAMETERS
###################################################
radius=float(argv[1])*DEG
try:qdeep=int(argv[2])
except:qdeep=False
try:qshow=int(argv[3])
except:qshow=False

print "Radius:",radius*RAD
fdata="scratch/points-r%.2e.data"%(radius*RAD)
fdataf="scratch/points-r%.2e-filtered.data"%(radius*RAD)

###################################################
#LOAD DATA
###################################################
data=np.loadtxt(fdata)
ss=data[:,2:]*RAD
ps=data[:,:2]

###################################################
#STATISTICS OF DISTANCES
###################################################
#"""
print "Statistics of distances..."
distmins=[]
distmaxs=[]
i=0
bad=0
print "Initial points: ",ss.shape[0]
ssgood=[]
psgood=[]
#timeIt()
dx=radius/PI
for s in ss:
    if (i%100)==0:
        print "Testing distance to %d..."%i
        #timeIt()

    l,b=s
    p=ps[i]

    ssearch=ss[np.abs(np.abs(ps[:,0])-np.abs(p[0]))<5*dx]
    #ssearch=ss

    i+=1

    dists=np.array([arcDistance(s*DEG,t*DEG) for t in ssearch])
    dists=dists[dists>1e-5]
    if min(dists)<=radius:
        bad+=1
    else:
        ssgood+=[s]
        psgood+=[p]
        distmins+=[min(dists)]
        distmaxs+=[max(dists)]

ssgood=np.array(ssgood)
psgood=np.array(psgood)
data=np.hstack((psgood,ssgood*DEG))
np.savetxt(fdataf,data)

print "Bad points: ",bad
print "Final points: ",ssgood.shape[0]

distmins=np.array(distmins)*RAD
distmaxs=np.array(distmaxs)*RAD

fig=plt.figure()
ax=fig.gca()
ax.hist(distmins)
fig.savefig("scratch/distance-grid-distrib.png")

fig=plt.figure()
ax=fig.gca()
ax.hist(distmaxs)
fig.savefig("scratch/distance-grid-distrib-max.png")

#"""

"""
data=np.loadtxt(fdataf)
ssgood=data[:,2:]
psgood=data[:,:2]
#"""

###################################################
#MAP
###################################################
print "Map of points..."
plt.close("all")
proj='robin'
map=drawMap(proj=proj,proj_opts=dict(lon_0=180),
            pars=[-45,0,45],mers=[0,45,90,135,180,225,270,315],
            pars_opts=dict(labels=[1,1,0,0],fontsize=8),
            mers_opts=dict(labels=[0,0,0,1],fontsize=8))
plotMap(map,np.mod(ssgood[:,0],360),ssgood[:,1],lw=0,
        marker='o',color='r',ms=2,mec='none')
plt.savefig("scratch/random-grid-blue-map.png")

###################################################
#3D POINTS MAP
###################################################
print "3D Map of points..."
ls=np.mod(ssgood[:,0],360.0)
bs=ssgood[:,1]

xs=np.cos(bs*DEG)*np.cos(ls*DEG)
ys=np.cos(bs*DEG)*np.sin(ls*DEG)
zs=np.sin(bs*DEG)

plt.close("all")
fig=plt.figure()
ax3d=plt3d(fig)
ax3d.set_aspect('equal')

u,v=np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x=np.cos(u)*np.sin(v)
y=np.sin(u)*np.sin(v)
z=np.cos(v)

ax3d.plot_wireframe(x, y, z, color="k")
ax3d.plot(xs,ys,zs,'o')
ax3d.plot_surface(x,y,z,rstride=1,cstride=1,color='c',alpha=1,linewidth=0)

plt.savefig("scratch/random-grid-blue-3d.png")
if qshow:plt.show()
