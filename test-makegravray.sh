#PBS -S /bin/bash
#PBS -j oe
#PBS -o log/makeagravray.log
cd $PBS_O_WORKDIR
python makeagravray.py "02/15/2013 3:20:34 UTC" deg locations.dat.temp rad directions.dat.temp 0 velocities.dat.temp "All temp" 80000
