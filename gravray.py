#############################################################
#EXTERNAL LIBRARIES
#############################################################
#REGULAR
from commands import getoutput as System
from sys import argv,exit
import numpy as np
from time import time

#GRAPHICAL
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D as plt3d,art3d

import visual as v

#DATABASE
import MySQLdb as mdb

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

TAB="\t"

###################################################
#DATABASE CONNECTION
###################################################
DATABASE="MinorBodies"
USER="minorbodies"
PASSWORD="347940"
CON=mdb.connect("localhost",USER,PASSWORD,DATABASE)
DB=CON.cursor()

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
    
ETIME=time()
TIME=time()
def timeIt():
    global TIME,ETIME
    ETIME=time()-ETIME
    total=time()-TIME
    print "[Elapsed: last %.2f msec, total %.3f]"%(1e3*ETIME,total)
    ETIME=time()

class dict2obj(object):
    def __init__(self,dic={}):self.__dict__.update(dic)
    def __add__(self,other):
        for attr in other.__dict__.keys():
            exec("self.%s=other.%s"%(attr,attr))
        return self

def mysqlSelect(selection="*",table="Bodies",condition="limit 100"):

    #FIELDS
    

    #QUERY
    sql="select %s from %s %s"%(selection,table,condition)
    DB.execute(sql)
    results=DB.fetchall()

    if selection=="*":
        #RECOVER ALL COLUMNS IN TABLE
        DB.execute("show columns from %s"%table)
        fields=DB.fetchall()
    else:
        #ONLY COLUMNS IN SELECTION
        fields=[(field,) for field in selection.split(",")]

    dresults=[]
    for result in results:
        row=dict()
        for i in xrange(len(fields)):
            field=fields[i][0]
            row[field]=result[i]
        dresults+=[row]

    return dresults
    
def listdict2matrix(listdict,keys=[]):
    n=len(listdict)
    m=len(listdict[0].keys())
    matrix=np.zeros((n,m))
    i=0
    for obj in listdict:
        if len(keys)==0:keys=obj.keys()
        matrix[i,:]=np.array([float(obj[key]) for key in keys])
        i+=1
    return matrix
    
    
