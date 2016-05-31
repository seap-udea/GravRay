def comparePairs(p1,p2):
    m=int(min((p2[0]-p1[0])*10000,(p2[1]-p1[1])*10000))
    return m

#SORT POINTS ACCORDING TO LOWEST COORDINATE
mincoord=np.array([min(x,y) for x,y in data])
isort=mincoord.argsort()

sdata1=np.array(list(reversed(sorted(data,cmp=comparePairs))))
sdata2=data[isort]

for p in sdata:
    print "Cartesian:",p
    print "Spherical:",car2sph(p)
exit(0)

