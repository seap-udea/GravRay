#!/bin/bash
nprocs=4
cmd="python generatedirections.py" 
#cmd="bash generatedirections.test" 
dirdata="util/data/directions"
basefile="direction-r*"

nrep=5
np=1
for angle in 30 20 15 12 10 7 5 3 1.5
do
    for i in $(seq 1 $nrep)
    do
	echo "Launching command $np..."
	command="$cmd $angle 1 0 sample$i"
	echo -e "\t$command"
	time $command >> log/prog-$np.log 2>&1 &
	((np++))
	if [ $np -eq $nprocs ];then
	    echo -e "\tWaiting..."
	    wait
	    np=0
	fi
    done
done
