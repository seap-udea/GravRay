#Low resolution
#python makeagravray.py "01/21/2019 04:41:37 UTC" deg locations.dat.moon rad util/data/directions-r3.00e+01.data 0 util/data/velocities.moon.regular50.dat "L21-J" 500

#High resolution
python makeagravray.py "01/21/2019 04:41:38 UTC" deg locations.dat.moon rad util/data/directions-r5.00e+00.data 0 util/data/velocities.moon.regular100.dat "L21-J" 500

#Test with a single integration
#python makeagravray.py "01/21/2019 04:41:37 UTC" deg locations.dat.moon deg directions.dat.moon 0 velocities.dat.moon "L21-J" 500
