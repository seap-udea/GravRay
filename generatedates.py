#!/usr/bin/env python
from gravray import *
import spiceypy as sp
sp.furnsh("kernels.txt")

#############################################################
#USAGE
#############################################################
usage="""Generate a list of dates from a repetition patern

   python <initial_date> <units> [<discrete> <delta1>,<delta2>,...]|[<grid> <span>,<delta>]|<linspace> <span>,<number>] <dates_file>

where:

   <initial_date>: Initial date.
                   Use format MM/DD/CCYY HH:MM:SS.dcm UTC-L.
   
   <units>: SEC, MIN, HOUR, DAY, YEAR

   Type of span: discrete, grid

     If 'discrete': 

         <delta1>,<delta2>,... : Times from last date.

     If 'grid':

        <span>: How much time from initial date.

        <delta>: Time among dates

     If 'linspace':

        <span>: How much time from initial date.

        <delta>: Number of dates
 
   <dates_file>: File to store list of dates
"""

#############################################################
#INPUT
#############################################################
iarg=1
try:
    inidate=argv[iarg];iarg+=1

    unitstr=argv[iarg];iarg+=1

    typespan=argv[iarg];iarg+=1
    spanpars=argv[iarg];iarg+=1

    datesfile=argv[iarg];iarg+=1
except:
    print usage
    exit(1)
  
#############################################################
#TIMES
#############################################################
et=sp.str2et(inidate)
unit=eval(unitstr)

if typespan=='discrete':
    pass
elif typespan=='grid':
    span,delta=(unit*float(f) for f in spanpars.split(","))
    dets=np.arange(0,span+delta,delta)
    times=et+dets
elif typespan=='linspace':
    unit=1
    span,num=(float(f) for f in spanpars.split(","))
    span*=unit
    dets=np.linspace(0,span,num+1)
    times=et+dets

#############################################################
#TIMES
#############################################################
for 
