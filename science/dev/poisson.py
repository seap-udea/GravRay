#-*- coding:utf-8 -*-
from gravray import *
from copy import copy

rand=np.random.rand
R=1
DEG=np.pi/180
RAD=1/DEG

def hav(theta):
    h=np.sin(theta/2)**2
    return h;

def arcDistance(p1,p2):
    dl=p2["l"]-p1["l"]
    db=p2["b"]-p1["b"]
    h=hav(db*DEG)+np.cos(p1["b"]*DEG)*np.cos(p2["b"]*DEG)*hav(dl*DEG)
    delta=2*np.arcsin(np.sqrt(h))*RAD
    return delta

def findClosest(candidate,points):
    for p in points:
        delta=arcDistance(candidate,p)
        if delta<candidate["delta"]:
            candidate["delta"]=delta

verbose=0
radius=2
points=[]

nest=int(2/(1-np.cos(radius*DEG)))
print "Estimated number of points:",nest

coords=[]
npoints=nest/3
for i in xrange(npoints):
    if (i%1)==0:print "Generating %d/%d..."%(i,npoints)

    best=None
    deltamax=0

    k=0
    maxk=10
    while True:

        if verbose:print "Points:",points
        if verbose:print "Best:",best

        candidate=dict(l=360*rand(),
                       b=90-np.arccos(1-2*rand())*180/np.pi,
                       delta=1E100)
        if verbose:print "Candidate:",candidate

        findClosest(candidate,points)
        if verbose:print "Delta minimum:",candidate["delta"]

        if best is None:
            if candidate["delta"]>2*radius:
                points+=[candidate]
        
        if candidate["delta"]>deltamax and candidate["delta"]>=2*radius and candidate["delta"]<1E100:
            best=copy(candidate)
            deltamax=candidate["delta"]

        if verbose:raw_input()
        k+=1
        if k>maxk:
            if best is None:
                maxk+=1
            else:
                print "\t%d darts used..."%k
                break
    
    if verbose:print "Deltamax: ",deltamax
    points+=[best]
    coords+=[[best["l"],best["b"]]]
    if verbose:raw_input

coords=np.array(coords)
np.savetxt("points.data",coords)

fig=plt.figure()
ax=fig.gca()

for p in points:
    ax.plot([p["l"]],[p["b"]],'ko')

ax.set_xlim((0,360))
ax.set_ylim((-90,90))
fig.savefig("scratch/my-bluenoise.png")

exit(0)

def poisson(k):
    delta=10
    points=[]
    geometries=[]
    best=dict(l=0,b=0,delta=0)
    for i in xrange(k):
        candidate=dict(l=360*rand(),
                       b=90-np.arccos(1-2*rand())*180/np.pi,
                       delta=1E100)
        findClosest(candidate,points)
        if best.delta==0 or candidate.radius>best.radius:
            best=copy(candidate)

        points+=[best]
        
poisson(10)

"""
// Poisson-disc sampling. Best of k candidates.
// Based on Mike Bostock’s http://bl.ocks.org/mbostock/1893974
function poisson(k) {
  var radius = 10,
      points = [],
      geometries = [],
      findClosest = finder(points, radius * 2);
  return function() {
    var best = null;

    // Create k candidates, picking the best (furthest away).
    for (var i = 0; i < k; ++i) {
      var candidate = {x: λ(Math.random()), y: 180 * Math.acos(Math.random() * 2 - 1) / Math.PI - 90};
      findClosest(candidate);
      if (!best || candidate.radius > best.radius) best = candidate;
    }

    best.radius = 5;
    points.push(best);
    geometries.push({type: "Point", coordinates: [best.x, best.y], id: nextId()});
    if (geometries.length > n) geometries.shift(), points.shift();
    return geometries;
  };
}

// Find the closest circle to the candidate.
function finder(points) {
  var arc = d3.geo.greatArc().target(Object);

  return function(candidate) {
    candidate.radius = Infinity;
    arc.source([candidate.x, candidate.y]);

    points.forEach(function(point) {
      var radius = Math.max(0, arc.distance([point.x, point.y]) * 180 / Math.PI - point.radius);
      if (radius < candidate.radius) candidate.radius = radius;
    });
  };
}
"""
