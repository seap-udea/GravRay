Gravitational Raytracing
========================

## A software by Jorge I. Zuluaga (C) 2016

Presentation
------------

*GravRay* is a package that implements the "Gravitational Raytracing"
method [1]. The method was originally devised by Jorge I. Zuluaga of
the Solar, Earth and Planetary Physics Group of the University of

The method was originally intended for studying the spatial
distribution of meteoroid and asteroid impacts on the Earth.  However
it can be adapted to other general purposes (impacts on other Solar
System bodies - the Moon, Mars, Jupiter, Saturn's Rings, etc., minor
bodies entrance to earth sphere of influence).

If you use or modify this package please cite:

>  [1] Zuluaga, J.I. and Sucerquia, M., Towards a theoretical
   determination of the geogpraphical probability distibution of
   meteoroid impacts, Planetary and Space Science, under submission
   (2016).

Get the package
---------------

To get the package directly from the github repository use:

```
   git clone http://github.com/seap-udea/GravRay.git
```

If you are a contributor and your key is already registered in GitHub
use:

```
   git clone git@github.com:seap-udea/GravRay.git
```

The package uses ~300 Mb of disk space (represented mainly by SPICE
Kernels which are used to calculate accurately the position of Solar
System objects).

Basic configuration
-------------------

Once cloned you need to perform the following basic actions (just once):

1. Create a copy of the template file ``compiler.in.temp`` as the file
   ``compile.in``:

   ```
       cp compile.in.temp compile.in
   ```

2. Set your system architecture.  Edit ``compiler`` file and uncomment
   your architecture (32 or 64 bits).  GravRay used precompiled
   versions of two of its main dependencies, namely CSPICE and GSL.

3. "Unpack" large kernels.  Large CSPICE kernels (those having more
   than 100 Mb) need to be splitted for distribution with git.  Before
   using the package "unsplit" them using:

   ```
       make unpack
   ```

4. Compile the codes.  For this you will need "g++`` installed in the
   machine.

   ```
       make
   ```

   If everything compiles as expected you're ready to use the package.

The MPC Database
----------------

Before using some modules of the package will also need a copy of the
Minor Planet Center database.

We have independently prepared a MySQL easy to install version of the
database that can be downloaded also from our github repository.

For downloading and installing the database go to the ``util`` directory
and run:

```
   git clone http://github.com/seap-udea/MPCDatabase.git
```

Then go to the MPCDatabase directory and run:

```
   cat database.sql user.sql | mysql -u root -p
```

You will need the password of root user in your mysql server.

Once the database is installed you need to populate it with the data.
We have automatized this procedure in a makefile:

```
   make restore
```

Structure of the package
------------------------

The package has a set of programs having specific functionalities.

Programs has easy to recall names such as ``whattimeisit`` or
``wherewillitbe``.

Below I summarize the functionality of the most important programs:

- ``whattimeisit``: convert from standard date format to ephemeris time
   and Julian Date.

- ``whereami``: compute the cartesian coordinates and velocity of an
   observer in the surface of the Earth with respect to the Solar
   System Barycenter in the ecliptic reference frame of the epoch.

- ``whereisit``: compute the position of a given major solar system
   object with respect to the SSB in the EclipticJ2000 reference
   frame.

- ``whereisthisasteroid``: the same as ``whereisit`` but for major
   Asteroids.

- ``whereisinsky``: calculates the position in the sky (RA, DEC, Az,
   Alt) of a given object with respect to a given observer in the
   surface of the Earth (precision of the results below 1 arcsec).

- ``wherewillitbe``: integrates the trajectory of a test particle in the
   scenario defined by the objects listed in ``objects.hpp``.

- ``scenarioof``: calculates the position of major solar system objects
   at the times defined in a ray trajectory (the output of
   wherewillitbe).

- ``scenarioof``: calculates the position of major solar system objects at
   the times defined in a ray trajectory.

- ``whereonearth``: calculates matrices defining certain regions of
   interest on the earth (angle with respect to velocity vector - apex
   = 1, antiapex = -1, cosine of ecliptic colatitude - ecliptic = 0,
   among others).

- ``throwrays``: this program will throw many rays starting in a given place
   in the world.  

The main routines and information required for all the programs are in
the master file ``gravray.cpp``.  The list of objects used on any
integration is defined in the ``objects.hpp`` file.  When you change
this file you should recompile the programs using ``make all``.  It will
only recompile the files depending on the objects file.

Testing the package
-------------------

Once the package is installed and configured you may test it by
running the programs with example inputs.  Here there are several
example tests you may run:

- Get the ephemeris time of a given date:

  ```
    ./whattimeisit.exe <date> <output_system>
  ```

  ```
    ./whattimeisit.exe "07/19/2015 00:00:00.000 UTC-5" ET
  ```

- Calculate the position of an observer on the Earth:

  ```
    ./whereami.exe <lat(deg)> <lon(deg)> <alt(m)> <elev(deg)> <azim(deg)> <vel(km/s)> <date>
  ```

  ```
    ./whereami.exe 6.2 -75.34 1450.0 45.0 0.0 1.0 "07/19/2015 00:00:00.000 UTC-5"
  ```

- Calculate the position of a major solar system object:

  ```
    ./whereisit.exe <object> [<date> | ET <time>]
  ```

  ```
    ./whereisit.exe MARS_BARYCENTER "07/19/2015 00:00:00.000 UTC-5"
  ```

  or

  ```
    ./whereisit.exe MARS_BARYCENTER ET 4.905360682e+08
  ```

- Calculate the position of an Asteroid:

  ```
    ./whereisit.exe <asteroid> [<date> | ET <time>]
  ```

  ```
    ./whereisthisasteroid.exe EROS "07/19/2015 00:00:00.000 UTC-5"
  ```

  or

  ```
    ./whereisthisasteroid.exe EROS ET 4.905360682e+08
  ```

  A list of the available asteroids is at the end of
  ``whereisthisasteroid.cpp`` source code file.

- Calculate position in the sky of a major solar system object as seen
  from a given place on the Earth:

  ```
    ./whereisinsky.exe <object> <lat(deg)> <lon(deg)> <alt(m)> <date>
  ```

  ```
    ./whereisinsky.exe MARS_BARYCENTER 6.2 -75.34 1450.0 "07/19/2015 00:00:00.000 UTC-5"
  ```

- Propagate the orbit of a body in a given gravitational scenario
  defined by the ``objects.hpp`` file:

  ```
    ./wherewillitbe.exe <time(ET seconds)> <x> <y> <z> <vx> <vy> <vz> <time(years)> <timesteps>
  ```

  ```
    ./wherewillitbe.exe 4.905360000e+08 -7.47594519221052825e+07 1.52138798910335392e+08 4.49404456025594100e+06 -2.69676500003677440e+01 -1.39288966833683254e+01 -5.76432883102505045e+00 +20.0 300
  ```

  This piece of code will propagate for instance the orbit of asteroid
  Eros 20 years in the future starting with its state vector of a
  given date.  The output will be stored in the file ``ray.dat``.  300
  data points will be considered for this calculation.

- Calculate the position and velocity of ``objects.hpp`` at the specific
  times defined in ``ray.dat``:

  ```
    ./scenarioof.exe <file>
  ```

  ```
    ./scenarioof.exe ray.dat
  ```

- Get the ecliptic colatitude and the colatitude of every point on
  Earth with respect to Earth's velocity at a given date:

  ```
    ./whereonearth.exe [ET|NOMAT] <date>|<ephemeris.time>
  ```

  ```
    ./whereonearth.exe "02/15/2013 03:20:34 UTC"
  ```

- Calculate the final orbital properties of a set of rays launched
  from a given location on Earth.  This is probably the more important
  program in the package since it is the basis for the GRT method.

  ```
    ./throwrays.exe <TDB> <lat> <lon> <alt> <input-file> <output-file>
  ```

  ```
    ./throwrays.exe 4.141704340e+08 54.456093 63.492323 8.234765e+04 initial.dat.temp elements.dat
  ```

  All input and outputs are given in km and km/s.

Stdout and Stderr
-----------------

You may control which information to display when using the programs
by redirecting the standard output and standard error.

Standard output shows a userfriendly output that can be redirectly to
/dev/null when you need plain information.

Standard error, on the other hand, handles the plain output of the
programas that can be manipulated for instance for wrapping scripts.

Thus you may run a program in two ways:

1. If you want the user friendly output:

   ```
   ./whattimeisit.exe "07/19/2015 00:00:00.000 UTC-5" TDB 2> /dev/null

   Input date: 07/19/2015 00:00:00.000 UTC
   ET = 4.905360682e+08
   Julian Date at ET = 2457222.500789
   TDB = 4.905360000e+08
   Julian Date at TDB = 2457222.500000
   Custom system (TAI) = 4.90536036000000000e+08
   ```

2. If instead you want only the plain information:

   ```
   ./whattimeisit.exe "07/19/2015 00:00:00.000 UTC-5" TDB > /dev/null

   ET,JD,DT,TDB,JDB,TAI
   4.905360682e+08
   2457222.500789
   68.18
   4.905360000e+08
   2457222.500000
   4.90536036000000000e+08
   ```

   You should notice that the plain information declares which
   information will be displayed below, in this case the time information
   returned by the program.

Testing the accuracy
--------------------

You can test the accuracy of the programs with different methods.

1. Comparing the results with NASA Horizons system.  For that purpose
   you may just run the ``tests/test-positions.sh`` script.

   ```
       bash tests/test-positions.sh "07/19/2015 00:00:00.000 UTC"
   ```

   This will create an output file ``scratch/test-positions.log`` with
   the position of the major planets at a given date.  Go to NASA
   Horizon system and compare your results with this system.

2. Integrate the orbit of an already existing object and then
   comparing the results with its expected position according to
   SPICE.  You can run this test with:

   ```
       python tests/test-integrator.sh Planet MOON
   ```

   or:

   ```
       python tests/test-integrator.sh Asteroid EROS
   ```

   This script will integrate the orbit of the given object using as
   force field the list of objects in ``objects.hpp``.  It is important
   that before running this script in the case of major planets, you
   disable in ``objects.hpp`` the respective planet. Not modification of
   this file is needed when is the case of an Asteroid.

   The test will produce the following files: ``dray.dat`` (difference
   between the integrated and the SPICE trajectory),
   ``dray-pos(vel).png``, ``dray-dist-pos(vel).png`` (different among
   the integrated position and velocities and the SPICE ones), and
   ``dray3d.png``, which is the trajectory of the ray in 3D.

   All these file are generated in the ``scratch`` directory.

3) Integrate forward then backward and orbit and calculate the
   difference. This is a classical test of an integrator.  You may run
   this test similarly as the previous tests:

   ```
       python tests/test-direction.py Planet MOON
   ```

   or:

   ```
       python tests/test-direction.py Asteroid EROS
   ```

   The script will generate two plots: ``dray-pos-direction.png`` and
   ``dray-vel-direction.png`` showing the relative error of position
   and velocities calculated in one direction then in the other.

Launching an asteroid
---------------------

Let's check that the package work with a particular impact that has
already happen on the Earth.  We will work with the Chelyabinsk impactor.

Determine the ecliptic coordinates of the impactor site:

```
./whereami.exe 54.456093 63.492323 8.234765e+04 17.664618 104.975030 2.045864e+01 "02/15/2013 3:20:34 UTC"
```

The output will look like:

```
Julian Date = 2456338.639282
Position observer w.r.t. ITRF93: 1.67989336300074933e+03 3.36821456733346031e+03 5.23342384710596252e+03
        Velocity of rotation w.r.t. ITRF93: -2.44943475473354683e-01 1.22165292778180709e-01 -9.01055435659073763e-18
        Velocity observer w.r.t. LOCAL (outwards): -5.03721069882287953e+00 -1.88319509387253916e+01 6.20806595218705137e+00
        Velocity observer w.r.t. LOCAL (inwards): 5.03721069882287953e+00 1.88319509387253916e+01 -6.20806595218705137e+00
        Velocity of motion w.r.t. ITRF93: 1.34092246260097987e+01 -1.53083252848880225e+01 -2.09852866334257726e+00
        Velocity total w.r.t. ITRF93: 1.31642811505364445e+01 -1.51861599921098414e+01 -2.09852866334257726e+00
Position observer w.r.t. J2000: -7.44564498708154247e+02 -1.30185961701703718e+03 6.26948335119827880e+03
        Velocity observer w.r.t. J2000: -1.66473939631271186e+01 9.50839416670916826e+00 -6.38563157989301633e+00
Position earth w.r.t. SSB J2000: -1.23266643210348308e+08 8.13492147654450387e+07 -9.73501068802284317e+03
        Velocity earth w.r.t. SSB J2000: -1.69443468835489668e+01 -2.49412909463430168e+01 1.48676424729513295e-03
Position observer w.r.t. SSB J2000: -1.23267387774847016e+08 8.13479129058280289e+07 -3.46552733682456437e+03
        Velocity observer w.r.t. SSB J2000: -3.35917408466760889e+01 -1.54328967796338485e+01 -6.38414481564572078e+00
--PLAIN--
TDB,JD,DT,ITRF93(6),ECJ2000(6)
4.141704340e+08
2456338.639282
67.19
+1.67989336300074933e+03 +3.36821456733346031e+03 +5.23342384710596252e+03 +1.31642811505364445e+01 -1.51861599921098414e+01 -2.09852866334257726e+00
-1.23267387774847016e+08 +8.13479129058280289e+07 -3.46552733682456437e+03 -3.35917408466760889e+01 -1.54328967796338485e+01 -6.38414481564572078e+00
```

Get the initial time in TDB and the coordinates of the impactor site
and propagate into the past for two years:

```
./wherewillitbe.exe 4.141704340e+08 -1.23267387774847016e+08 +8.13479129058280289e+07 -3.46552733682456437e+03 -3.35917408466760889e+01 -1.54328967796338485e+01 -6.38414481564572078e+00 -2.0 100
```

The output will be available in the file ``ray.dat``.  In particular
you should be interested in the last line:

```
tail -n 1 ray.dat

+3.51055234000000000e+08  -2.78806305352541864e+08  -1.86545885718426019e+08  -3.16928385508019067e+07  -3.95340660787797837e-01  -1.71759717456054979e+01  -1.48971663716463243e+00  +3.36952161883870482e+08  +1.72449863816401887e+01  +1.05023695277718261e+08  +6.12003711316022270e-01  +5.84849118490910058e+00  +3.26515953696134829e+02  +1.02986331263678224e+02  +8.14434356167327422e+01  +3.51055234000000000e+08  +1.32692605000000000e+11
```

Whose columns are: #1:t, 2:x, 3:y, 4:z, 5:vx, 6:vy, 7:vz, 8:r, 9:v,
10:q, 11:e, 12:i, 13:W, 14:w, 15:M, 16:t0, 17:mu.

Python Scripts and GRT Analysis
-------------------------------

In order to perform a Ray Tracing analysis of a given site or a given
geographical region we have prepared a set of python scripts that
using the aofrementioned programs performed all the required dynamical
analysis.

Below we describe the scripts:

- ``generatevelocities.py``: generate a set of impact velocities
  following a given distribution.

  ```
  Usage:

        python generatevelocities.py [<nsample>] [<source>] [<method>]

  Where:

	nsample: number of velocities to be generated.

	source: which source of objects velocities do you want to use.
        	Available sources are: velimp (impact velocities),
        	velinf (velocity at infinite), velast (velocities with
        	respect to SSB)

	method: generation method, either random (following PDF) or
        	regular (spaced proportionally to probability)

  Example:

	python generatevelocities.py 10 velimp regular
  ```

- ``generatedirections.py``: This script generate random directions
   obeying a blue noise distribution, i.e. the distance between the
   points is larger than a given radius.  

  ``` 
  Usage:

     python generatedirections.py <radius in degrees> [<perform a deep analysis?>]

  Where:

	radius: minimum distance in degrees between generated points.

	qdeep: 1 if you want that the final points generated be
               analysed deeply to see if some of them are separated by
               distances lower than the minimum.  All points (not only
               the closer ones) will be analysed.  If qdeep=0 the
               vicinity analysis will be done only over the closer
               points.

  The method was devised by Jorge I. Zuluaga and could be very
  inefficient for very small values of the radius (<3 degrees).

  Example:

    python generatevelocities.py 10.0 1 1

  ```

  This directions could be used either for initial geographical
  positions or for initial directions in the sky.

- ``generateinitialconditions.py``: This script generate random directions
   obeying the random noise distribution, i.e. the distance between
   the points is larger than a given radius.

  ```
  Usage:

      python generateinitialconditions.py <file.dirs> <file.vels> <file.initial>

  Where:

	<file.dirs>: file with random directions to be used as
		     Azimuths (longitude) and latitudes.  It can be
		     generated using generatedirections.py
		     script. There are a set of precalculated
		     directions that could be used here and that are
		     in the util/data directory

	<file.vels>: file with the list of velocities that will be
	             used in the simulation.

	<file.initial>: output file cointaining the initial conditions
  	 		of the simulation (a matrix with Az, h and v).
  
  ```
   

- 

- ``throwaray.py``: throw a particle from a given location on Earth and
   at a given velocity.

  ```
  Usage:

	python throwaray.py <lat> <lon> <alt> <elev> <azim> <vimp> <date> <tspan> [npoints]

  Where:

    <lat>,<lon>: latitude and longitude in degrees.
    <alt>: altitude in meters.
    <elev>,<azim>: elevation and azimuth in degrees.
    <vimp>: impact velocity in km/s
    <date>: date in format MM/DD/CCYY HH:MM:SS.dcm UTC-L
    <tspan>: integration time (in years)
    <npoints>: sampling points

  Example: This example launch to space the Chelyabinsk impactor.  The
  negative velocity is chosen because the object is impacting;
  correspondingly the total time is also chosen negative.

    python throwaray.py 54.456093 63.492323 8.234765e+04 17.664618 104.975030 -2.045864e+01 "02/15/2013 3:20:34 UTC" -2.0

  ```
