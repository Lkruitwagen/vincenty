Vincenty
========

Calculate direct and inverse Vincenty methods. The direct method calculates the destination coordinates given an origin coordinate, a geodesic (great circle) distance, and a departure azimuth angle. The inverse method calculates a geodesic distance given origin and destination coordinates. The direct method also returns the arrival azimuth angle, and the inverse method returns both the departure and arrival azimuth angles. 

Vincenty's method solves for geodesic distances on an ellipsoid and is widely used in geographic information systems. It is much more accurate than spherical-earth assumptions and is accurate within 1mm. The implementation given here uses the WGS84 ellipsoid, a standard reference geometry for Earth.

Example: distance between Boston and New York City
--------------------------------------------------

.. code:: python

   >>> from vincenty import *
   >>> boston = (42.3541165, -71.0693514)
   >>> newyork = (40.7791472, -73.9680804)
   >>> V_inv(boston, newyork)
   #distance [km], azimuth_dept [deg], azimuth_arr [deg]
   (298.396057, 235.0838926194711, 233.1602005520544)
   >>> V_inv(boston, newyork, miles=True)
   #distance [miles], azimuth_dept [deg], azimuth_arr [deg]
   (185.414713, 235.0838926194711, 233.1602005520544)
   >>> V_inv(newyork, boston)
   #distance [km], azimuth_dept [deg], azimuth_arr [deg]
   (298.396057, 53.1602005520544, 55.0838926194711)
   >>> V_dir(boston, 298.396057, 235.0838926194711)
   #arrival_point (lat [deg],lon [deg]), azimuth_arr [deg]
   ((40.779147202568396, -73.96808039548996), 233.16020055500005)


References
----------

* https://en.wikipedia.org/wiki/Vincenty's_formulae
* https://en.wikipedia.org/wiki/World_Geodetic_System
