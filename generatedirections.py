from gravray import *

###################################################
#PARAMETERS
###################################################
radius=float(argv[1])*DEG
print "Radius:",radius*RAD
fdata="scratch/points-r%.2e.data"%(radius*RAD)

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

    #SKIP DOUBLE OCCUPIED CELLS
    if ng>=2:continue

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
        if not qaccept:break

    #============================== 
    #CHECK DISTANCE TO POINTS 
    #============================== 
    if qaccept:
        if (n%100)==0:
            print "Point %d/%d/rej.%d accepted in cell "%(n,0.7*fullfill,nrej),ix,iy,":",p
        ngrid[ix,iy]+=1
        pgrid[ix][iy]+=[p]
        nrej=0
        n+=1

    #============================== 
    #FILL COEFFICIENT
    #============================== 
    cond=ngrid>0
    fill=int(ngrid[cond].sum())
    if n>0.7*fullfill:
        print "End by completion"
        break
    if nrej>0.5*fullfill:
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
            l,b=s
            ix,iy=xy2ij(p[0],p[1])
            ncells=[(ix,iy)]+genCells(l,b)
            qaccept=True
            for ic,jc in ncells:
                cell=pgrid[ic][jc]
                if len(cell)==0:continue
                for t in cell:
                    dist=normedArcDistance(p,t)
                    if dist<=radius and dist>0:
                        qaccept=False
                        break
                if not qaccept:break
            if qaccept:
                ps+=[p]
                ss+=[s]

ps=np.array(ps)
ss=np.array(ss)
data=np.hstack((ps,ss))
print "Final number:",data.shape[0]
np.savetxt(fdata,data)
