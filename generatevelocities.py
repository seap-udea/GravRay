from gravray import *

###################################################
#INPUTS
###################################################
nsample=int(argv[1])
velfile="scratch/velocities-%.4e.data"%(nsample)

###################################################
#LOAD REQUIRED INFORMATION
###################################################
vp=np.loadtxt("util/velocity-pdf.data")
vmin=vp[:,0].min();vmax=vp[:,0].max();

###################################################
#GENERATE VELOCITIES
###################################################
#RANDOM 
#vs=generateVelocities("util/velocity-cum.data",nsample)

#UNIFORMLY DISTRIBUTED IN CUMMULATIVE 
du=1./nsample
u=np.linspace(du,1,nsample)
velcum=np.loadtxt("util/velocity-cum.data")
ifvelcum=interp1d(velcum[:,1],velcum[:,0])
vs=ifvelcum(u)

###################################################
#SAVE VELOCITIES
###################################################
print "Saving sample with %d impact velocities..."%nsample
np.savetxt(velfile,vs)

###################################################
#HISTOGRAM
###################################################
nbins=vp.shape[0]
h,x=np.histogram(vs,nbins,(vmin,vmax))
fig=plt.figure()
ax=fig.gca()
ax.plot(x[:-1],h,'+')
ax.plot(vs,np.ones_like(vs),'ko',ms=4)
ax.plot(vp[:,0],vp[:,1]*nsample,'-')
fig.savefig("scratch/velocities-pdf-sample.png")
