.. _inputfiles_:

*****************
Input File Format
*****************

.. moduleauthor:: Matthew Burger <mburger@stsci.edu>

Input files are plain text files in the form: ::

    category.parameter = setting

Lines in the input file that can not be parsed in this manner are ignored.
Comments can be entered with either a "," or a "#". Everything in a line
after a comment character is ignored. There are currently seven categories
that can be set: geometry, surface_interaction,
forces, spatial_dist, speed_dist, angular_dist, and options. The required
parameters for each category are not fixed; i.e., which paramters are needed
depends somewhat on the settings chosen. Below, all possible parameters for
each category are defined. Input files are case insensitive.

Geometry
========

The geometry can be defined either with a timestamp (i.e., determine the
geometry values at a defined epoch), or without a timestamp (i.e., by
specifying important values). If running the model from MESSENGERdata.model(),
the geometry is determined from from the data and geometry settings in the
input file are ignored.

Geometry With Time Stamp
------------------------

geometry.planet [Required]
    Central planet for the model. This must be an object that orbits the Sun.

geometry.StartPoint [Optional]
    Object from which packets are ejected. This must be an object in the
    planetary system (the planet or one of its moons).
    Default = geometry.planet.

geometry.objects [Optional]
    Objects to include in calculations given as comma-separated list of
    bodies in the planetary system. For example, if
    `geometry.objects = Jupiter, Io`, the gravity effects of the other moons
    would not be included, nor would collisions with their surfaces. Default
    is to include all defined objects for the system.

geometry.startttime [Required]
    Starting time for the model run in ISOT format (YYYY-MM-DDTHH:MM:SS).
    The true anomaly angle, subsolar point, and orbital position of moons
    are determined using `SPICE <https://naif.jpl.nasa.gov/naif/toolkit.html>`_.

Geometry Without Time Stamp
---------------------------

geometry.planet [Required]
    Central planet for the model. This must be an object that orbits the Sun.

geometry.StartPoint [Optional]
    Object from which packets are ejected. This must be an object in the
    planetary system (the planet or one of its moons).
    Default = geometry.planet.

geometry.objects [Optional]
    Objects to include in calculations given as comma-separated list of
    bodies in the planetary system. For example, if
    `geometry.objects = Jupiter, Io`, the gravity effects of the other moons
    would not be included, nor would collisions with their surfaces. Default
    is to include all defined objects for the system.

geometry.phi [Required if planet has moons.]
    Orbital phase of each included moon relative to the Sun in radians given
    as a comma-separated list.
    Measured from 0 rad to 2π rad where 0 rad is superior conjunction and
    π/2 rad is over the planet's dawn terminator. The number of values must be
    equal to the number of non-planet objects included.

geometry.subsolarpoint [Optional]
    The sub-solar longitude and latitude over the planet's surface in radians
    given as comma-separated values. For Jupiter, use System-III central
    meridian longitude. Sub-solar latitude isn't used for anything currently,
    but could in the future be used to include effects of the planetary system's
    tilt relative to the Sun.

geometry.TAA [Required]
    Planet's True Anomaly Angle in radians. This is used to determine the
    planet's distance and radial velocity relative to the Sun.

SurfaceInteraction
==================

The SurfaceInteraction class defines interactions between packets and body
surfaces. The available parameters depend on the interactions desired.
If no values are provided, 100% sticking is assumed.

Constant Sticking Coefficient
-----------------------------

surfaceinteraction.stickcoef [Optional]
    Sticking coefficient to be used uniformly across the body's surface.
    For complete surface sticking, set `surfaceinteraction.stickcoef = 1.`.
    For no sticking (100% of packets are reemitted from the surface, set
    `surfaceinteraction.stickcoef = 0`. Default = 1.

surfaceinteraction.accomfactor [Required if stickcoef < 1]
    Surface accommodation factor. 1 = Fully accommodated to local surface
    temperature. 0 = Elastic reemission.

Temperature Dependent Sticking Coefficient
------------------------------------------

The sticking coefficient follows the functional form (Yakshinskiy & Madey 2005):

.. math::
    S(T) = A_0 e^{A_1 T} + A_2

where the coefficients are species dependent. For Na,
:math:`A_0=1.57014, A_1=-0.006262, A_2=0.1614157`.

surfaceinteraction.sticktype [Required]
    Set `surfaceinteraction.sticktype = temperature dependent`.

surfaceinteractions.accomfactor [Required]
    Surface accommodation factor. 1 = Fully accommodated to local surface
    temperature. 0 = Elastic reemission.

surfaceinteractions.A [Optional]
    Comma separated values for the coeffeicients.
    Default = 1.57014, 0.006262, 0.1614157. (Ideally the defaults will be
    species dependent, but I only have values for Na.)



Sticking Coefficient from a Surface Map
---------------------------------------

surfaceinteraction.type [Required]
    Set `surfaceinteraction.sticktype = from map`.

surfaceinteraction.sticking_mapfile [Required]
    Path to the file containing a map of the sticking coeficient. The format
    for the map has not been determined.

surfaceinteractions.accomfactor [Required]
    Surface accommodation factor. 1 = Fully accommodated to local surface
    temperature. 0 = Elastic reemission.

Forces
======

The Forces class determines which forces are included in the simulation.
Currently, the model only includes gravity and radiation pressure. If
no forces are set in the input file both are included by default.

forces.gravity [Optional]
    True to include gravity; False to exclude. Default = True.

forces.radpres [Optional]
    True to include radiation pressure; False to exclude. Default = True

SpatialDist
===========

The SpatialDist class specifies the initial spatial distribution of packets
in the system. Currently, three spatial distribution types are defined, all of
which place packets over the surface (or exobase) of *geometry.StartingPoint*.
More distributions may defined upon request.

**Coordinate Systems**

The coordinate system used for the object's latitude and longitude depends
on whether the packets are ejected from a planet or a moon. For planets, a
solar-fixed coordinate system is used where the longitude increases in the
positive direction from the sub-solar point (noon) point to dusk point: ::

    sub-solar (noon) point = 0 rad = 0°
    dusk point = π/2 rad = 90°
    anti-solar (midnight) point = π rad = 180°
    dawn point = 3π/2 rad = 270°

For satellites, the coordinate system is planet-fixed from the sub-planet
point increasing positive through the leading point: ::

    sub-planet point = 0 rad = 0°
    leading point = π/2 rad = 90°
    anti-planet point = π rad = 180°
    trailing point = 3π/2 rad = 270°

Latitude ranges from -π/2 rad to π/2 rad for the south pole to the north pole.
All angular values are given in radians in the input file.

Uniform Surface
---------------

Distribute packets randomly across a region of the surface or exobase with
a uniform probability distribution.

spatialdist.type [Required]
    Set `spatialdist.type = uniform`.

spatialdist.longitude [Optional]
    Longitude range on the surface to place packets in radians given as
    *long0, long1* where :math:`0 \leq long0,long1 \leq 2\pi`. If *long0* >
    *long1*, the region wraps around. Default = 0, 2π.

spatialdist.latitude [Optional]
    Latitude range on the surface to place packets in radians given as
    *lat0, lat1* where :math:`-\pi/2 \leq lat0 \leq lat1 \leq \pi/2`.

spatialdist.exobase [Optional]
    Location of the exobase in units of the starting point's radius.
    Default = 1.

To eject all packets from a single point, set *long0 = long1* and
*lat0 = lat1*; i.e., to eject all packets from the sub-solar point of a planet,
set: ::

    spatialdist.longitude = 3.14159,3.14159
    spatialdist.latitude = 0,0

Spatial Distribution from a Surface Map
---------------------------------------

Distribute packets according to a probability distribution given by a
pre-defined surface map.

spatialdist.type [Required]
    Set `spatialdist.type = surface map`.

spatialdist.mapfile [Required]
    Path to the file containing the surface map.

spatialdist.coordsystem [Optional]
    Coordinate system for the map. Options are solar-fixed, planet-fixed,
    or moon-fixed. Default = solar-fixed for planets and planet-fixed for
    moons. For a surface composition map that rotates with the planet, use
    planet-fixed.

spatialdist.exobase [Optional]
    Location of the exobase in units of the starting point's radius.
    Default = 1.

Surface-Spot Spatial Distribution
---------------------------------

Distribute packets with a spatial distribution that drops off exponentially
from a central point.

spatialdist.type [Required]
    Set `spatialdist.type = surface spot`.

spatialdist.longitude [Required]
    Longitude of the source center in radians.

spatialdist.latitude [Required]
    Latitude of the soruce center in radians.

spatialdist.sigma [Required]
    Angular e-folding width of the source in radians.

spatialdist.exobase [Optional]
    Location of the exobase in units of the starting point's radius.
    Default = 1.

SpeedDist
=========

The SpeedDist class defines the one-dimensional initial speed distribution
of the packets. Currently implemented speed distributions are gaussian,
Maxwellian, sputtering, and flat. More can be added upon request.

Gaussian (Normal) distribution
------------------------------

Packets speeds are chosen from a normal distribution. See
`numpy.random.normal
<https://docs.scipy.org/doc/numpy-1.16.0/reference/generated/numpy.random.normal.html#numpy.random.normal>`_
for more information on the implementation.

speeddist.type [Required]
    Set `speeddist.type = gaussian`

speeddist.vprob [Required]
    Mean speed of the distribution in km/s.

speeddist.sigma [Required]
    Standard deviation of the distribution in km/s.

Maxwellian Distribution
-----------------------

Packet speeds are chosen from a Maxwellian distribution given by:

.. math::
    :nowrap:

    \begin{eqnarray*}
    f(v) & \propto & v^3 \exp(-v^2/v_{th}^2) \\
    v_{th}^2 & = & 2Tk_B/m
    \end{eqnarray*}

speeddist.type [Required]
    Set `speeddist.type = maxwellian`

speeddist.temperature [Required]
    Temperature of the distribution in K. Set `speeddist.temperature = 0` to
    use a pre-defined surface temperature map (Not implemented yet).

Sputtering Distribution
-----------------------

Packet speeds are chosen from a sputtering distribution in the form:

.. math::
    :nowrap:

    \begin{eqnarray*}
    f(v) & \propto & \frac{v^{2\beta + 1}}{(v^2 + v_b^2)^\alpha} \\
    v_b & = & \left(\frac{2U}{m} \right)^{1/2}
    \end{eqnarray*}

speeddist.type [Required]
    Set `speeddist.type = sputtering`

speeddist.alpha [Required]
    :math:`\alpha` parameter.

speeddist.beta [Required]
    :math:`\beta` parameter.

speeddist.U [Required]
    Surface binding energy in eV.

Flat Distribution
-----------------

Packet speeds are uniformly distributed between *vprob - delv/2* and
*vrpob + delv/2*. Setting `speeddist.delv = 0` gives a monoenergetic
distribution.

speeddist.type [Required]
    Set `speeddist.type = flat`

speeddist.vprob [Required]
    Mean speed of the distribution in km/s.

speeddist.delv [Required]
    Full width of the distribution in km/s.

AngularDist
===========

The AngularDist class defines the initial angular distribution of packets.
The options are radial and isotropic. More distributions can be added upon
request. If not given, an isotropic distribution into the outward facing
hemisphere is assumed.

Radial Distribution
-------------------

Packets are ejected radially from the surface.

angulardist.type [Required]
    Set `angulardist.type = radial`.

Isotropic Distribution
----------------------

Packets are ejected isotropically into the outward facing hemisphere (if the
packets are starting from the surface) or the full hemisphere.
`angulardist.type` is not given, an isotropic distribution is assumed and
all other options are ignored (i.e., altitude and azimuth can not be specified).

angulardist.type [Optional]
    Set `angulardist.type = isotropic`.

angulardist.altitude [Optional]
    Used to limit the altitude range of the distribution. Given as a
    comma-separated list of *altmin, altmax* in radians measured from the
    surface tangent to the surface normal.

angulardist.azimuth [Optional]
    Used to limit the azimuth range of the distribution. Given as a
    comma-separated list of *az0, az1* in radians. This should be measured with
    azimuth = 0 rad pointing to north, but I'm not sure if it actually works.
    Use of this option is not recommended.

Options
-------

The Options class sets runtime options that don't fit into other categories.

options.endtime [Required]
    The total simulated runtime for the model. Generally chosen to be several
    times the lifetime of the species.

options.species [Required]
    The species to be simulated.

options.lifetime [Optional]
    The lifetime due to ionization or dissociation of the species in seconds.
    If `options.lifetime = 0`, the lifetime is computed based on available
    ionization and dissociation reactions. If `options.lifetime > 0`, the
    lifetime is constant throughout the system. If `options.lifetime < 0`,
    the lifetime is assumed to be the photo-lifetime and no loss occurs in
    the geometric shadow. Default = 0 (use available reactions).

options.outer_edge [Optional]
    Distance from *geometry.startpoint* to simulate in object radii. Default =
    infinite; i.e., no outer edge is given to the simulation.

options.step_size [Optional]
    Time step size for the simulation in seconds. Set `options.step_size = 0`
    for variable step size. Default = 0 (variable step size). If step_size is
    non-zero, the number of steps to be run is endtime/step_size + 1.

options.resolution [Optional]
    Relative precision of the simulation. Default = :math:`10^{-4}`.
    This is ignored if *options.step_size* is set.
