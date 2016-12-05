#!/bin/bash
points="scratch/points-r1.00e+01-filtered.data"

lons=$(awk '{printf "%s\n",$3*180/3.1415926538}' $points)
for lon in $lons
do
    echo $lon
done
