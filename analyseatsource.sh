#!/bin/bash
grtdir=$1
locals="$grtdir/locals.dat"
tini=$(date +%s)
for ray in $grtdir/rays-*.data
do
    python analyseatsource.py $locals $ray
done
tend=$(date +%s)
time=$((tend-tini))
print "Total execution time:$time secs"
