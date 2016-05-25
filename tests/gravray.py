from commands import getoutput as System
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D as plt3d
from sys import argv
import numpy as np

#############################################################
#MACROS
#############################################################
norm=np.linalg.norm

#############################################################
#CONSTANTS
#############################################################
#TIME UNITS
MIN=60.0
HOUR=60.0*MIN
DAY=24*HOUR
YEAR=365.25*DAY

#UNITCONVERSION
UL=1.49597870691E8 #km
UV=2.97824659599673905e+01 #UL/UT

#FACTOR CONVERSION
STATECONV=np.array([UL,UL,UL,UV,UV,UV])

#OBJECTS
OBJECTS={"0":"Sun",
         "1":"Mercury",
         "2":"Venus",
         "3":"Earth",
         "4":"Moon",
         "5":"Mars",
         "6":"Jupiter",
         "7":"Saturn",
         "8":"Uranus",
         "9":"Neptune"}

#############################################################
#ROUTINES
#############################################################
def vec2str(vec,frm="%+.17e "):
    string=""
    for comp in vec:
        string+=frm%comp
    return string
    
def out2state(out,ini=0,end=-1):
    state=np.array([float(x) for x in out.split(" ")[ini:end]])
    return state

def getScenario(sfile):
    scenario=np.loadtxt(sfile)
    ncols=scenario.shape[1]
    nobjs=(ncols-1)/7
    obj=dict()
    obj["ts"]=scenario[:,0]
    for i in xrange(nobjs):
        stride=1+7*i
        dataobj=scenario[:,stride:stride+7]
        iobj=int(dataobj[0,0])
        obj["%d"%iobj]=dataobj[:,1:]

    return obj

def out2dict(output):
    parts=output.split("\n")
    comps=parts[0].split(",")
    i=0
    output=dict()
    for part in parts[1:]:
        if '(' in comps[i]:
            values=part.split(" ")
            output[comps[i]]=np.array([float(value) for value in values[:-1]])
        else:
            output[comps[i]]=float(part)
        i+=1
    return output
    
