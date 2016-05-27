from gravray import *

if int(argv[1])>0:
    props="Perihelion_dist,e,i"
    listdict=mysqlSelect(props,
                         condition="where NEO_flag and e>0 and e<1")
    elements=listdict2matrix(listdict,keys=props.split(","))

def plotDistribution(el):
    qes=el[:,0];qmin=qes.min();qmax=qes.max()
    ees=el[:,1];emin=ees.min();emax=ees.max()
    ies=np.log10(el[:,2]);imin=ies.min();imax=ies.max()
    aes=qes/(1-ees)

    #INTERVALS
    qlow=0.6;qup=1.3
    elow=0.2;eup=0.7
    ilow=0.0;iup=2.0
    
    #CUTS
    iconst=1.2;di=0.1*10 # Logarithmic
    qconst=1.2;dq=0.1*100
    econst=0.5;de=0.05*100

    fig=plt.figure(figsize=(12,8))
    fig3d=plt.figure(figsize=(8,8))

    ax3d=fig.add_subplot(221,projection='3d')
    axae=fig.add_subplot(222)
    axai=fig.add_subplot(223)
    axei=fig.add_subplot(224)

    # CUTS
    condi=np.abs(ies-iconst)<di
    q=[qmin,qmin,qmax,qmax,qmin]
    e=[emin,emax,emax,emin,emin]
    i=[iconst]*5
    verts=[zip(q,e,i)]
    ax3d.add_collection3d(art3d.Line3DCollection(verts,color='black'))

    condq=np.abs(qes-qconst)<dq
    q=[qconst]*5
    e=[emin,emin,emax,emax,emin]
    i=[imin,imax,imax,imin,imin]
    verts=[zip(q,e,i)]
    ax3d.add_collection3d(art3d.Line3DCollection(verts,color='black'))

    conde=np.abs(ees-econst)<de
    q=[qmin,qmin,qmax,qmax,qmin]
    e=[econst]*5
    i=[imin,imax,imax,imin,imin]
    verts=[zip(q,e,i)]
    ax3d.add_collection3d(art3d.Line3DCollection(verts,color='black'))

    ax3d.plot(qes,ees,ies,'o',mec='none',ms=1)


    scplot=dict(color='k',mec='none',ms=2,alpha=0.0)

    H,xe,ye=np.histogram2d(qes[condi],ees[condi],bins=50,normed=True)
    scale=(qup-qlow)/(eup-elow)
    img=axae.imshow(H.transpose(),origin='lower',interpolation='bicubic',extent=(qmin,qmax,emin,emax),aspect=scale)
    axae.plot(qes[condi],ees[condi],'o',**scplot)
    axae.text(0.05,0.90,"i = (%.2f,%.2f)"%(10**(iconst-di),10**(iconst+di)),transform=axae.transAxes)

    H,xe,ye=np.histogram2d(qes[conde],ies[conde],bins=50,normed=True)
    scale=(qup-qlow)/(iup-ilow)
    img=axai.imshow(H.transpose(),origin='lower',interpolation='bicubic',extent=(qmin,qmax,imin,imax),aspect=scale)
    axai.plot(qes[conde],ies[conde],'o',**scplot)
    axai.text(0.05,0.90,"e = (%.2f,%.2f)"%(econst-de,econst+de),transform=axai.transAxes)

    H,xe,ye=np.histogram2d(ees[condq],ies[condq],bins=50,normed=True)
    scale=(eup-elow)/(iup-ilow)
    img=axei.imshow(H.transpose(),origin='lower',interpolation='bicubic',extent=(emin,emax,imin,imax),aspect=scale,cmap='spectral')
    axei.plot(ees[condq],ies[condq],'o',**scplot)
    axei.text(0.05,0.90,"q = (%.2f,%.2f)"%(qconst-dq,qconst+dq),transform=axei.transAxes)

    #DECORATION
    ax3d.set_xlabel("q")
    ax3d.set_ylabel("e")
    ax3d.set_zlabel("i")

    #DECORATION
    axae.set_xlabel("q")
    axae.set_ylabel("e")
    axae.set_xlim((qlow,qup))
    axae.set_ylim((elow,eup))

    axai.set_xlabel("q")
    axai.set_ylabel("log i")
    axai.set_xlim((qlow,qup))
    axai.set_ylim((ilow,iup))

    axei.set_xlabel("e")
    axei.set_ylabel("log i")
    axei.set_xlim((elow,eup))
    axei.set_ylim((ilow,iup))

    fig.savefig("scratch/distribution.png")

