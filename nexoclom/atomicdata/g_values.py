"""``g_values`` - Routines related to g-values and radiation pressure
The g-value is the the product of the solar flux at the dopler-shifted
emission wavelength and the scattering probability per atom. See
`Killen, R.M. et al., Icarus 209, 75–87, 2009.
<http://dx.doi.org/10.1016/j.icarus.2010.02.018.>`_ for details on calculating
g-values for important species in Mercury's atmosphere.

The radiation acceleration is given by
:math:`a_{rad} = h g/m \lambda`,
where h is Plank's constant, g is the g-value as a function of radial
velocity, m is the mass of the accelerated species, and λ is the wavelength
of the absorbed photon.
"""
import numpy as np
import pandas as pd
import astropy.units as u
from astropy import constants as const
from nexoclom.atomicdata.atomicmass import atomicmass
from nexoclom.utilities.database_connect import database_connect
from nexoclom.math import interpu


class gValue:
    r"""Class containing g-value vs. velocity for a specified atom and
    transition.

    **Parameters**

    sp
        atomic species

    wavelength
        Wavelength of the transition. Default=None.

    aplanet
        Distance from the Sun. Can be given as an astropy quantity with
        distance units or as a float assumed to be in AU. Default = 1 AU

    **Class Attributes**

    species
        The input species

    wavelength
        The input wavelength

    aplanet
        The input aplanet

    velocity
        Radial velocity deviation relative to the Sun in km/s.
        Positive values indicate
        motion away from the Sun. Given as a numpy array of astropy quantities

    g
        g-value as function of velocity in units 1/s.
    """
    def __init__(self, sp, wavelength=None, aplanet=1*u.au):
        self.species = sp
        if wavelength is None:
            assert 0
            # waves = pd.read_sql(
            #     f'''SELECT DISTINCT wavelength
            #         FROM gvalues
            #         WHERE species='{self.species}' ''', con)
            # self.wavelength = [w * u.AA for w in waves.wavelength]

        try:
            self.wavelength = wavelength.to(u.AA)
        except:
            self.wavelength = wavelength * u.AA

        try:
            self.aplanet = aplanet.to(u.au)
        except:
            self.aplanet = aplanet * u.au

        with database_connect() as con:
            gvalue = pd.read_sql(
                f'''SELECT refpt, velocity, g
                    FROM gvalues
                    WHERE species='{self.species}' and
                          wavelength='{self.wavelength.value}' ''', con)

        if len(gvalue) == 0:
            self.velocity = np.array([0., 1.])*u.km/u.s
            self.g = np.array([0., 0.])/u.s
            print(f'Warning: g-values not found for species = {sp}')
        elif len(gvalue) == 1:
            self.velocity = np.array(gvalue.velocity[0])*u.km/u.s
            self.g = (np.array(gvalue.g[0])/u.s *
                      gvalue.refpt[0]**2/self.aplanet.value**2)
        else:
            assert 0, 'Multiple rows found.'


class RadPresConst:
    r"""Class containing radial acceleration vs. velocity for a specified atom.

    **Parameters**

    sp
        atomic species

    aplanet
        Distance from the Sun. Can be given as an astropy quantity with
        distance units or as a float assumed to be in AU. Default = 1 AU

    database
        Database containing solar system information. Default =
        `thesolarsystem` which probably shouldn't be overridden.

    **Class Attributes**

    species
        The input species

    aplanet
        The input distance from the Sun

    velocity
        Radial velocity deviation relative to the Sun in km/s.
        Positive values indicate
        motion away from the Sun. Given as a numpy array of astropy quantities

    accel
        Radial acceleration vs. velocity with units km/s**2.
    """
    def __init__(self, sp, aplanet):
        self.species = sp
        try:
            self.aplanet = aplanet.value * u.au
        except:
            self.aplanet = aplanet * u.au

        # Open database connection
        with database_connect() as con:
            waves = pd.read_sql(f'''SELECT DISTINCT wavelength
                                    FROM gvalues
                                    WHERE species='{sp}' ''', con)

        if len(waves) == 0:
            self.v = np.array([0., 1.])*u.km/u.s
            self.accel = np.array([0., 0.])*u.km/u.s**2
            print(f'Warning: g-values not found for species = {sp}')
        else:
            self.wavelength = [w*u.AA for w in waves.wavelength]

            gvals = [gValue(sp, w, aplanet) for w in self.wavelength]

            # Complete velocity set
            allv = []
            for g in gvals:
                allv.extend(g.velocity.value)
            allv = np.unique(allv) * u.km/u.s

            # Interpolate gvalues to full velocity set and compute rad pres
            rpres = np.zeros_like(allv)/u.s
            for g in gvals:
                g2 = interpu(allv, g.velocity, g.g)
                rpres_ = const.h/atomicmass(sp)/g.wavelength * g2
                rpres += rpres_.to(u.km/u.s**2)
                
            self.velocity = allv
            self.accel = rpres
