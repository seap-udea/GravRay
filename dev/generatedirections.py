from gravray import *
rand=np.random.rand
n=0

def hav(theta):
    h=np.sin(theta/2)**2
    return h;

def arcDistance(p1,p2):
    dl=p2["l"]-p1["l"]
    db=p2["b"]-p1["b"]
    h=hav(db*DEG)+np.cos(p1["b"]*DEG)*np.cos(p2["b"]*DEG)*hav(dl*DEG)
    delta=2*np.arcsin(np.sqrt(h))*RAD
    return delta

def car2sph(p):
    l=np.mod(180*p[0],360)
    b=90.0-np.arccos(p[1])*RAD
    return l,b

def distCartesian(p1,p2):
    return norm(p1-p2)

def distSpherical(p1,p2):
    s=car2sph(p1)
    p1=dict(l=s[0],b=s[1])
    s=car2sph(p2)
    p2=dict(l=s[0],b=s[1])
    return arcDistance(p1,p2)

#############################################################
#GET DISTANCE DISTRIBUTION (SPHERICAL)
#############################################################
data=np.loadtxt("data_5e-2.data")
cradius=0.02*2
radius=cradius*RAD
#SORT BY COORDINATE
sdata=np.array(sorted(data,key=lambda p:(p[0],p[1])))
cond=1-np.abs(sdata[:,0])>cradius/100
sdata=sdata[cond]

print "Size original:",sdata.shape[0]
v=0
print "Maximum radius:",radius
i=0
distmins=[]
distcmins=[]
for point in sdata:
    if (i%100)==0:print "Point %d..."%i
    x=point[0]
    y=point[1]
    if v:print "Point:",point,car2sph(point)
    cond=(np.abs(sdata[:,0]-x)<2*cradius)*((np.abs(sdata[:,1]-y)<2*cradius))
    neighbors=sdata[cond]
    dists=np.array([distSpherical(point,neighbor) for neighbor in neighbors])
    distcs=np.array([distCartesian(point,neighbor) for neighbor in neighbors])
    try:
        mindist=min(dists[dists>0])
        mindistc=min(distcs[distcs>0])
        if v:print "Neighbors:",[car2sph(p) for p in neighbors]
        if v:print "Minimum distance:",mindist
        if mindist<radius:
            if v:print "Removing point:",sdata[i,:]
            sdata=np.delete(sdata,i,0)
        else:
            distmins+=[mindist]
            distcmins+=[mindistc]
    except:
        pass
    i+=1
    if v:raw_input()
        
distmins=np.array(distmins)
distcmins=np.array(distcmins)
print "Size final:",sdata.shape[0]

fig=plt.figure()
ax=fig.gca()
ax.hist(distmins)
fig.savefig("scratch/distance-distrib.png")

fig=plt.figure()
ax=fig.gca()
ax.hist(distcmins)
fig.savefig("scratch/distance-cartesian-distrib.png")

r1=[];r2=[]
for x,y in sdata:
    if x>0 and y>0 or 1:
        r1+=[x];r2+=[y]
r1=np.array(r1);r2=np.array(r2)
f=np.mod(np.pi*r1*RAD,360)
q=np.arccos(r2**(1/(n+1.)))*RAD
l=90-q

proj='robin'
map=drawMap(proj=proj,proj_opts=dict(lon_0=180),
            pars=[-45,0,45],mers=[0,45,90,135,180,225,270,315],
            pars_opts=dict(labels=[1,1,0,0],fontsize=8),
            mers_opts=dict(labels=[0,0,0,1],fontsize=8))
plotMap(map,f,l,lw=0,marker='o',color='r',ms=3,mec='none')
plt.savefig("scratch/random-directions-blue-map.png")

exit(0)

n=0

"""
cond=(np.abs(data[:,1])<1-0.1)
data=data[cond]
#"""
distmins=[]
for point in data:
    distmin=1e100
    for other in data:
        dist=distSpherical(point,other)
        #dist=distCartesian(point,other)
        if dist<distmin and dist>0:
            distmin=dist
    distmins+=[distmin]
distmins=np.array(distmins)

print distmins.mean()

fig=plt.figure()
ax=fig.gca()
ax.hist(distmins)
fig.savefig("scratch/distance-distrib.png")

r1=[];r2=[]
for x,y in data:
    if x>0 and y>0 or 1:
        r1+=[x];r2+=[y]
r1=np.array(r1);r2=np.array(r2)
f=np.mod(2*np.pi*(1-r1)/2*RAD,360)
q=np.arccos(r2**(1/(n+1.)))*RAD
l=90-q

proj='robin'
map=drawMap(proj=proj,proj_opts=dict(lon_0=180),
            pars=[-45,0,45],mers=[0,45,90,135,180,225,270,315],
            pars_opts=dict(labels=[1,1,0,0],fontsize=8),
            mers_opts=dict(labels=[0,0,0,1],fontsize=8))
plotMap(map,f,l,lw=0,marker='o',color='r',ms=3,mec='none')
plt.savefig("scratch/random-directions-blue-map.png")
exit(0)

#############################################################
#GET DISTANCE DISTRIBUTION (CARTESIAN)
#############################################################
data=np.loadtxt("data_5e-2.data")

distmins=[]
for point in data:
    distmin=1e100
    for other in data:
        dist=norm(point-other)
        if dist<distmin and dist>0:
            distmin=dist
    distmins+=[distmin]
distmins=np.array(distmins)

fig=plt.figure()
ax=fig.gca()
ax.hist(distmins)
fig.savefig("scratch/distance-distrib.png")
exit(0)

#############################################################
#DIRECTIONS OBTAINED WITH SAMPLES FROM POISSON DISKS
#############################################################
N=1000
n=0

data=np.loadtxt("data.data")
r1=[];r2=[]
for x,y in data:
    if x>0 and y>0 or 1:
        r1+=[x];r2+=[y]
r1=np.array(r1);r2=np.array(r2)

q=np.arccos(r1**(1/(n+1.)))*RAD
f=np.mod(np.pi*r2*RAD,360)
l=90-q

plt.close("all")
#proj='cyl'
#proj='ortho'
#proj='hammer'
proj='robin'
map=drawMap(proj=proj,proj_opts=dict(lon_0=180),
            pars=[-45,0,45],mers=[0,45,90,135,180,225,270,315],
            pars_opts=dict(labels=[1,1,0,0],fontsize=8),
            mers_opts=dict(labels=[0,0,0,1],fontsize=8))
plotMap(map,f,l,lw=0,marker='o',color='r',ms=3,mec='none')
plt.savefig("scratch/random-directions-blue-map.png")
exit(0)

#############################################################
#RANDOM DIRECTIONS
#############################################################

#Phong-like distribution 
#See: http://www.cs.rutgers.edu/~decarlo/readings/mcrt-sg03c.pdf, p. 23
#For a uniform distribution over a sphere use n=0
r1=rand(N)
r2=rand(N)

q=np.arccos((1-2*r1)**(1/(n+1.)))*RAD
f=2*np.pi*r2*RAD
l=90.-q

fig=plt.figure()
ax=fig.gca()

ax.plot(f,l,'o')

ax.set_xlim((0,360))
ax.set_ylim((-90,90))
fig.savefig("scratch/random-directions.png")

plt.close("all")
map=drawMap(proj='hammer',proj_opts=dict(lon_0=180),
            pars=[-45,0,45],mers=[0,45,90,135,180,225,270,315],
            pars_opts=dict(labels=[1,1,0,0],fontsize=8),
            mers_opts=dict(labels=[0,0,0,1],fontsize=8))
plotMap(map,f,l,lw=0,marker='o',color='r',ms=3,mec='none')
plt.savefig("scratch/random-directions-map.png")

#############################################################
#DIRECTIONS OBTAINED WITH SAMPLES FROM POISSON DISKS
#############################################################
data=np.loadtxt("data.data")
r1=data[:,0]
r2=data[:,1]

print r1.shape
exit(0)

r1=r1[(r1>0)*(r2>0)]
r2=r2[(r1>0)*(r2>0)]

q=np.arccos((1-2*r1)**(1/(n+1.)))*RAD
f=2*np.pi*r2*RAD
l=90.-q

plt.close("all")
map=drawMap(proj='hammer',proj_opts=dict(lon_0=180),
            pars=[-45,0,45],mers=[0,45,90,135,180,225,270,315],
            pars_opts=dict(labels=[1,1,0,0],fontsize=8),
            mers_opts=dict(labels=[0,0,0,1],fontsize=8))
plotMap(map,f,l,lw=0,marker='o',color='r',ms=3,mec='none')
plt.savefig("scratch/random-directions-blue-map.png")
