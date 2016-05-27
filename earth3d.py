from numpy import *
from visual import *
from sys import argv,exit

nfiles=len(argv)-1
if nfiles>0:
    files=argv[1:]
else:
    files=["impactor.dat"]
    nfiles=1

print "%d files provided: "%nfiles,files

def phy2vis(rp):
    rv=vector(rp[1],rp[2],rp[0])
    return rv

def maketex(fileimg,highres=True):
    if highres:
        try:
            import Image
            im=Image.open(fileimg)
            im.resize((128,128),Image.ANTIALIAS)
            materials.saveTGA("texture",im)
            tex=materials.texture(data=im,mapping="spherical")
        except:
            tex=materials.earth
    else:
        tex=materials.earth
    return tex

#PHYSICAL CONSTANTS
RE=6371

#"""
#SCENE
scene.center=(0,0,0)
scene.forward=(-1,-1,-1)
scene.width=800
scene.height=600
scene.ratio=scene.width/scene.height
scene.range=(2*RE,2*RE,2*RE)
scene.ambient=0.8
scene.lights=[distant_light(direction=(0,0,-1),color=color.gray(1.0))]

#AXES
AR=1.2*RE
curve(pos=[(0,0,0),(AR,0,0)],color=color.white)
curve(pos=[(0,0,0),(0,AR,0)],color=color.white)
curve(pos=[(0,0,0),(0,0,AR)],color=color.white)
label(text='y',pos=(AR,0,0),color=color.white,align='center',box=False,opacity=0)
label(text='z',pos=(0,AR,0),color=color.white,align='center',box=False,opacity=0)
label(text='x',pos=(0,0,AR),color=color.white,align='center',box=False,opacity=0)
#EARTH
texture=maketex("textures/Earth.png",highres=True)
earth=sphere(pos=(0,0,0),radius=RE,material=texture,opacity=1.0)
#"""

fname=argv[1]

mycolor=[color.white,color.blue,color.red,color.yellow]
ncolors=len(mycolor)
for i in xrange(0,nfiles):
    fname=files[i]
    print "Loading info from: ",fname
    data=loadtxt(fname)

    if i==0:
        #PLOT REFERENCE POINTS
        xre=data[1,:]
        rr=xre[0:3]
        rv=phy2vis(rr)
        theta=xre[3]
        points(pos=[rv],size=10)
        label(text='Geodetic Origin',pos=rv,color=color.white,align='center',box=False,opacity=0)
        earth.rotate(angle=theta*pi/180,axis=(0,1,0))
        
        xrp=data[2,:]
        rr=xrp[0:3]
        rv=phy2vis(rr)
        psi=xrp[3]
        points(pos=[rv],size=10)
        label(text='North Pole',pos=rv,color=color.white,align='center',box=False,opacity=0)
        earth.rotate(angle=-psi*pi/180,axis=(0,0,1))
        
    #IMPACT POINT AND VELOCITY
    xe=data[0,:]
    rp=xe[0:3]
    vp=xe[3:7]
    r1=rp
    r2=array(r1)+array(vp)/mag(vp)*RE/5
    
    #PLOT POINT
    rv=phy2vis(rp)
    points(pos=[rv],size=10,color=mycolor[i%ncolors])

    #PLOT VELOCITY
    rv1=phy2vis(r1)
    rv2=phy2vis(r2)
    curve(pos=[rv1,rv2],size=10,color=mycolor[i%ncolors])
