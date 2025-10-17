"""Determine distance and radial velocity relative to Sun."""
import numpy as np
from scipy.misc import derivative
import astropy.units as u
# from astroquery.jplhorizons import Horizons
from nexoclom.solarsystem.SSObject import SSObject


def planet_dist(planet_, taa=None, time=None):
    """ Given a planet and either a TAA or a time, return distance from and
    radial velocity relative to the Sun.
    """
    if isinstance(planet_, str):
        planet = SSObject(planet_)
        
        if planet.object is None:
            return None
    elif isinstance(planet_, SSObject):
        planet = planet_
    else:
        raise TypeError('solarsystemMB.planet_dist',
                        'Must give a SSObject or a object name.')

    if time is not None:
        raise NotImplementedError
        # et = spice.str2et(time.isot)
        # posvel, lt = spice.spkezr(planet.object, et, 'J2000',
        #                           'LT+S', 'Sun')
        #
        # position = np.array(posvel[0:3])*u.km
        # r = np.sqrt(np.sum(position**2))
        #
        # velocity = np.array(posvel[3:])*u.km/u.s
        # v_r = np.sum(position*velocity)/r
        # r = r.to(u.au)
    elif taa is not None:
        a = planet.a
        eps = planet.e

        # make sure taa is in radians. If not a quantity, assume it is.
        if isinstance(taa, type(1*u.s)):
            taa_ = taa.to(u.rad)
        elif type(taa) in (int, float, np.float64):
            taa_ = taa * u.rad
        else:
            raise TypeError('taa must be a number or angle quantity')

        if eps > 0:
            # determine r
            r = a * (1-eps**2)/(1+eps*np.cos(taa_))
            period = planet.orbperiod.to(u.s)
            
            # determine v_r = dr/dt
            time= np.linspace(0, 1, 1000)*period.value
            time= np.concatenate([np.array([time[0]-time[1]]), time])*u.s
            
            mean_anomaly = np.linspace(0, 2*np.pi, 1000)
            mean_anomaly = np.concatenate(
                [np.array([mean_anomaly[0]-mean_anomaly[1]]), mean_anomaly])*u.rad
            
            true_anomaly = (mean_anomaly +
                            (2*eps - eps**3/4)*np.sin(mean_anomaly)*u.rad +
                            5/4 * eps**2 * np.sin(2*mean_anomaly)*u.rad +
                            13/12 * eps**3 * np.sin(3*mean_anomaly)*u.rad)
            r_true = a * (1-eps**2)/(1+eps*np.cos(true_anomaly))
            drdt = (r_true[1:] - r_true[:-1])/(time[1:] - time[:-1])
            v_r = np.interp(taa_, true_anomaly[1:], drdt.to(u.km/u.s))
        else:
            r, v_r = a, 0.*u.km/u.s
    else:
        print('Neither a time nor a true anomaly was given.')
        return None

    return r, v_r
