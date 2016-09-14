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
        Ptot=data[:,7].sum()/(1.0*Ninitial)
        dataensamble+=[[lat,lon,Ptot]]
        print fname,lat,lon,Ptot
        break
        
    dataensamble=np.array(dataensamble)
    print dataensamble
    exit(0)

