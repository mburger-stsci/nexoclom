import numpy as np
from nexoclom.atomicdata import atomicmass
import astropy.units as u
import astropy.constants as const


def sputdist(velocity, U, alpha, beta, species):
    mspecies = atomicmass(species)
    v_b = np.sqrt(2*U/mspecies)
    v_b = v_b.to(u.km/u.s)
    f_v = velocity**(2*beta+1) / (velocity**2 + v_b**2)**alpha
    f_v /= np.max(f_v)
    return f_v.value


def MaxwellianDist(velocity, temperature, species):
    vth2 = 2*temperature*const.k_B/atomicmass(species)
    vth2 = vth2.to(u.km**2/u.s**2)
    f_v = velocity**3 * np.exp(-velocity**2/vth2)
    f_v /= np.max(f_v)
    return f_v.value
