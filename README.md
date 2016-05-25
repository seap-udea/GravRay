Gravitational Raytracing
========================

A software by Jorge I. Zuluaga (C) 2016

Presentation
------------

*GravRay* is a package that implements the "Gravitational Raytracing"
method. The method was originally devised by Jorge I. Zuluaga of the
Solar, Earth and Planetary Physics Group of the University of
Antioquia.  The method was intended for studying the spatal
distribution of meteoroid and asteroid "impacts" on the surface of the
Earth.  However it can be adapted to other general purposes.

If you use or modify this package please cite:

>   Zuluaga, J.I. and Sucerquia, M., Towards a theoretical
>   determination of the geogpraphical probability distibution of
>   meteoroid impacts, Space Science Journal, in press (2016).

Get the package
---------------

To get the package use:

```
   git clone http://github.com/seap-udea/GravRay.git
```

If you are a contributer and your key is already registered in Git Hub
use:

```
   git clone git@github.com:seap-udea/GravRay.git
```
  
The package uses ~300 Mb of disk space (represented mainly by SPICE
Kernels used to calculate accurately the position of Solar System
objects).

Basic configuration
-------------------

Once cloned you need to perform the following basic actions (just once):

1. Set your system architecture.  Edit "compiler.in" file and
   uncomment your architecture (32 or 64 bits).  GravRay used
   precompiled versions of two of its main dependencies, namely CSPICE
   and GSL.

2. "Unpack" large kernels.  Large CSPICE kernels (those having more
   than 100 Mb) need to be splitted for distribution with git.  Before
   using the package "unsplit" them using:

```
       make unpack
```

3. Compile the codes.  For this you will need "g++" installed in the
   machine.
   
```
       make
```

If everything compiles as expected you're ready to use the package.

Structure of the package
------------------------

The package has a set of programs having specific functionalities.
Programs has easy to recall names such as "whattimeisit" or
"wherewillitbe". 

Below I summarize the functionality of the most important programs:

- *whattimeisit*: convert from standard date format to ephemeris time
   and Julian Date.

- *whereami*: compute the cartesian coordinates and velocity of an
   observer in the surface of the Earth with respect to the Solar
   System Barycenter in the ecliptic reference frame of the epoch.

- *whereisit*: compute the position of a given major solar system
   object with respect to the SSB in the EclipticJ2000 reference
   frame.

- *whereisthisasteroid*: the same as "whereisit" but for major
   Asteroids.

- *wherewillitbe*: integrates the trajectory of a test particle in the
   scenario defined by the objects listed in *objects.hpp*.

- *whereisinsky*: calculates the position in the sky (RA, DEC, Az,
   Alt) of a given object with respect to a given observer in the
   surface of the Earth.

- *scenario*: calculates the position of major solar system objects at
   the times defined in a ray trajectory.

The main routines and information required for all the programs are in
the master file *gravray.cpp*.  The list of objects used on any
integration is defined in the *objects.hpp* file.  When you change
this file you should recompile the programs using "make all".  It will
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
2C

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

- Calculate position in the sky of a major solar system object as seen
  from a given place on the Earth:

```
    ./whereisinsky.exe <object> <lat(deg)> <lon(deg)> <alt(m)> <date>
```

```
    ./whereisinsky.exe MARS_BARYCENTER 6.2 -75.34 1450.0 "07/19/2015 00:00:00.000 UTC-5"
```

- Propagate the orbit of a body in a given gravitational scenario
  defined by the *objects.hpp* file:

```
    ./wherewillitbe.exe <time(ET seconds> <x> <y> <z> <vx> <vy> <vz> <time(years)> <timesteps>
```

```
    ./wherewillitbe.exe 4.905360000e+08 -7.47594519221052825e+07 1.52138798910335392e+08 4.49404456025594100e+06 -2.69676500003677440e+01 -1.39288966833683254e+01 -5.76432883102505045e+00 +20.0 300
```

  This piece of code will propagate for instance the orbit of asteroid
  Eros 20 years in the future starting with its state vector of a
  given date.  The output will be stored in the file *ray.dat*.  300
  data points will be considered for this calculation.

- Calculate the position and velocity of *objects.hpp* at the specific
  times defined in *ray.dat*:

```
    ./scenario.exe <file>
```

```
    ./scenario.exe ray.dat
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

2) If instead you want only the plain information:

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
   you may just run the "tests/test-positions.sh" script.
   
```
       bash tests/test-positions.sh
```

   This will create an output file "scratch/test-positions.log" with
   the position of the major planets at a given date.  Go to NASA
   Horizon system and compare your results with this system.

2. Integrating the orbit of an already existing object and then
   comparing the results with its expected position according to
   SPICE.  You can run this test with:

```
       bash tests/test-integrator.sh Planet MOON
```

   or:

```
       bash tests/test-integrator.sh Asteroid EROS
```

   This script will integrate the orbit of the given object using as
   force field the list of objects in *objects.hpp*.  It is important
   that before running this script in the case of major planets, you
   disable in *objects.hpp* the respective planet. Not modification
   is needede when is the case of an Asteroid.

   The test will produce the following files: dray.dat (difference
   between the integrated and the SPICE trajectory),
   dray-pos(vel).png, dray-dist-pos(vel).png (different among the
   integrated position and velocities and the SPICE ones), and
   dray3d.png, which is the trajectory of the ray in 3D.

3) Integrating forward then backward and orbit and calculate the
   difference. This is a classical test of an integrator.  You may run
   this test similarly as the previous tests:

```
       bash tests/test-direction.sh Planet MOON
```

   or:
   
```
       bash tests/test-direction.sh Asteroid EROS
```

   The script will generate two plots: *dray-pos-direction.png* and
   *dray-vel-direction.png* showing the relative error of position and
   velocities calculated in one direction then in the other.

Launching an asteroid
---------------------

Let's check that the package work with a particular impact that has
already happen on the Earth.  We will work with the Chelyabinsk impactor.

1. Determine the ecliptic coordinates of the impactor site.

```
    ./whereami.exe 54.456093 63.492323 8.234765e+04 17.664618 104.975030 2.045864e+01 "02/15/2013 3:20:34 UTC"
```

2. Get the initial coordinates of the impactor