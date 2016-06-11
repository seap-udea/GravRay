from gravray import *

###################################################
#USAGE
###################################################
usage="""Generate random directions over a sphere following a blue noise
distribution (points are separated by a distance larger than a given
minimum).

python generatedirections.py <radius in degrees> [<perform a deep analysis?>]

Where:

  radius: minimum distance in degrees between generated points.

  qdeep: 1 if you want that the final points generated be analysed
         deeply to see if some of them are separated by distances
         lower than the minimum.  All points (not only the closer
         ones) will be analysed.  If qdeep=0 the vicinity analysis
         will be done only over the closer points.

The method was devised by Jorge I. Zuluaga and could be very
inefficient for very small values of the radius."""

###################################################
#INPUT
###################################################
iarg=1
try:
    radius=float(argv[iarg])*DEG
    iarg+=1
except:
    print usage
    exit(1)

try:
    qdeep=int(argv[iarg])
    iarg+=1
except:qdeep=False

print "Radius:",radius*RAD

fdata="scratch/directions-r%.2e-unfiletered.dat"%(radius*RAD)
fdataf="directions-r%.2e.dat"%(radius*RAD)

###################################################
#PARAMETERS
###################################################
#Maximum number of points accepted in a cell
MAXOCCUPY=5

#

###################################################
#CREATE GRID IN THE UNIT-SQUARE
###################################################

#========================================
#X-AXIS DIVISION
#========================================
dx=2*radius/(180*DEG)
xs=np.arange(-1,1+dx,dx)
xs=xs[xs<=1]
nx=xs.shape[0]

#========================================
#Y-AXIS DIVISION
#========================================
q=0;y=0
ys=[y]
while q<=90*DEG:
    q+=2*radius
    y=np.sin(q)
    ys+=[y]
ys[-1]=1.0
ys=np.array(ys)
ys=np.concatenate((-ys[::-1],ys[1:]))
ny=ys.shape[0]

ies=np.arange(nx-1)
js=np.arange(ny-1)
delta=np.array([normedArcDistance([xs[i],ys[j]],[xs[i+1],ys[j+1]]) for i,j in zip(ies,js)])*RAD

#========================================
#GRID OF DARTS
#========================================
ngrid=np.zeros((nx-1,ny-1))
pgrid=[[[] for j in xrange(ny-1)] for i in xrange(nx-1)]
fullfill=2*(nx-1)*(ny-1)

#========================================
#GRID ROUTINES
#========================================
def xy2ij(x,y):
    """
    Get the (i,j) of the nearest point in the grid to x,y

    Example: xy2ij(-0.5,0.5)
    """
    dx=x-xs[:-1]
    dy=y-ys[:-1]
    ix=ies[dx>=0][-1]
    iy=js[dy>=0][-1]
    return ix,iy

def ij2next(l,b,di,dj):
    """
    Get the (i,j) of a point at (l+di*r,b+dj*r)

    Example: ij2next(175.0,86.0,+1,+1)
    """
    lp=l+2*di*radius
    bp=b+2*dj*radius

    if bp>=90*DEG:
        bp=180*DEG-bp
        lp=180*DEG+lp
    if bp<=-90*DEG:
        bp=-180*DEG-bp
        lp=lp-180*DEG

    lp=np.mod(lp,360*DEG)
    if lp>180*DEG:lp=lp-360*DEG

    xp,yp=sph2car([lp,bp])
    ip,jp=xy2ij(xp,yp)
    
    return ip,jp

def genCells(l,b):
    cells=[
        ij2next(l,b,+0,-1),
        ij2next(l,b,+0,+1),
        ij2next(l,b,+1,-1),
        ij2next(l,b,+1,+0),
        ij2next(l,b,+1,+1),
        ij2next(l,b,-1,-1),
        ij2next(l,b,-1,+0),
        ij2next(l,b,-1,+1)
    ]
    return cells

###################################################
#THROW DARTS
###################################################
v=0
imax=len(ies)-1
jmax=len(js)-1

n=1
nrej=0
np.random.seed(1)
while True:
    #============================== 
    #INITIALIZE CELLS
    #============================== 
    iesr=ies
    jsr=js
    
    #============================== 
    #RANDOM POSITION
    #============================== 
    p=1-2*rand(2)
    x,y=p
    
    #CONVERT TO SPHERICAL
    l,b=car2sph([x,y])

    #============================== 
    #GRID CELL
    #============================== 
    ix,iy=xy2ij(x,y)
    ng=ngrid[ix,iy]

    #SKIP POINT IF CELL IS OCCUPIED BY MAXOCCUPY POINTS
    if ng>=MAXOCCUPY:continue

    #============================== 
    #NEIGHBOR CELLS
    #============================== 
    ncells=[(ix,iy)]+genCells(l,b)

    #============================== 
    #CHECK DISTANCE TO POINTS 
    #============================== 
    qaccept=True
    for ic,jc in ncells:
        cell=pgrid[ic][jc]
        if len(cell)==0:continue
        for t in cell:
            dist=normedArcDistance(p,t)
            if dist<=radius:
                qaccept=False
                nrej+=1
                break
        if not qaccept:break

    #============================== 
    #CHECK DISTANCE TO POINTS 
    #============================== 
    if qaccept:
        if (n%100)==0:
            print TAB,"Point %d/%d/rej.%d accepted in cell "%(n,
                                                              fullfill,
                                                              nrej),ix,iy,":",p
        ngrid[ix,iy]+=1
        pgrid[ix][iy]+=[p]
        nrej=0
        n+=1

    #============================== 
    #FILL COEFFICIENT
    #============================== 
    cond=ngrid>0
    fill=int(ngrid[cond].sum())
    if n>MAXOCCUPY*fullfill:
        print "End by completion"
        break
    if nrej>fullfill:
        print "End by rejection"
        break

print "%d points accepted"%n
print "Fill fraction: %d/%d"%(n,fullfill)
print "Last rejection rate: %d"%(nrej)

###################################################
#SAVE POINTS
###################################################
ps=[]
ss=[]
for i in xrange(nx-1):
    for j in xrange(ny-1):
        for p in pgrid[i][j]:
            s=car2sph(p)
            ps+=[p]
            ss+=[s]

ps=np.array(ps)
ss=np.array(ss)
data=np.hstack((ps,ss))
print "Final number:",data.shape[0]
np.savetxt(fdata,data)

###################################################
#FILTER POINTS
###################################################
print "Filter points..."
print "Initial points: ",ss.shape[0]

ssgood=[]
psgood=[]
distmins=[]
distmaxs=[]
dx=radius/PI
i=0
bad=0

sscomp=ss
indesa=np.arange(sscomp.shape[0])
for s in ss:
    if (i%100)==0:
        print TAB,"Testing distance of %d..."%(i)

    l,b=s
    p=ps[i]

    if not qdeep:
        cond=np.abs(np.abs(ps[:,0])-np.abs(p[0]))<=5*dx
        ssearch=sscomp[cond]
        indes=indesa[cond]
    else:
        indes=indesa
        ssearch=sscomp
    i+=1

    dists=np.array([arcDistance(s,t) for t in ssearch])
    isort=dists.argsort()

    #Check if the point is too close to other points
    if dists[isort[1]]<=radius:
        bad+=1
        ipresent=indes[isort[0]]
        iclosest=indes[isort[1]]
        """
        print "Distances:",sorted(dists*RAD)[:5]
        print "Sorted indexes: ",isort[:5]
        print "Closest point: ",iclosest
        print "Bad point at:",s*RAD
        print "Bad point at (by index):",sscomp[ipresent]*RAD
        print "Closest point:",sscomp[iclosest]*RAD
        print "Distance recalculated:",arcDistance(s,sscomp[iclosest])*RAD
        """
        sscomp=np.delete(sscomp,ipresent,0)
        indesa=np.arange(sscomp.shape[0])
    else:
        ssgood+=[s]
        psgood+=[p]
        distmins+=[min(dists[dists>1e-5])]
        distmaxs+=[max(dists[dists>1e-5])]

ssgood=np.array(ssgood)
psgood=np.array(psgood)
data=np.hstack((psgood,ssgood*RAD))
np.savetxt(fdataf,data)

print "Bad points: ",bad
print "Final points: ",ssgood.shape[0]

###################################################
#PLOT INFORMATION ABOUT POINTS
###################################################

#========================================
#DISTANCE STATISTICS
#========================================
distmins=np.array(distmins)*RAD
distmaxs=np.array(distmaxs)*RAD

fig=plt.figure()
ax=fig.gca()
ax.hist(distmins)
fig.savefig("scratch/distmin-distrib.png")

fig=plt.figure()
ax=fig.gca()
ax.hist(distmaxs)
fig.savefig("scratch/distmax-distrib.png")

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

plt.savefig("scratch/directions-3d.png")
if qshow:plt.show()
