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
