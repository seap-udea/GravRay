#High resolution
#python makeagravray.py "02/01/2019 18:17:00 UTC" deg locations.dat.cuba rad util/data/directions-r5.00e+00.data 0 util/data/velocities.earth.regular100.dat "Cuba" 30000

#Single
#python makeagravray.py "02/01/2019 18:17:00 UTC" deg locations.dat.cuba deg directions.dat.cuba 0 velocities.dat.cuba "Cuba_single" 4000

#Error
python makeagravray.py "02/01/2019 18:17:00 UTC" deg locations.dat.cuba_error deg directions.dat.cuba_error 0 velocities.dat.cuba_error "Cuba_error" 4000
