#!/bin/bash
NP=$1;shift
NR=$(cat gravray-include.sh | wc -l)
NSPLIT=$((NR/NP))
echo "Splitting $NR in $NP chunks of $NSPLIT rays..."
cd tmp/
rm gravray*
split ../gravray-include.sh -l $NSPLIT -d gravray
cd -
rm orbits*.dat
i=0
for file in tmp/gravray*
do
    echo "Running process $i..."
    bash $file $i 2>&1 |tee gravray$i.log > /dev/null &
    ((i++))
done
wait

(echo '#1:t                      2:x                       3:y                       4:z                       5:vx                      6:vy                      7:vz                      8:r                       9:v                       10:q                      11:e                      12:i                      13:W                      14:w                      15:M                      16:t0                     17:mu';cat orbits[0-9]*.dat) > orbits.dat
