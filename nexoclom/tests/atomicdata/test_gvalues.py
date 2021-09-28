import os.path
from IPython.terminal.embed import embed
import numpy as np
import pandas as pd
import astropy.units as u
from astropy import constants as const
from nexoclom.atomicdata import atomicmass, gValue, RadPresConst
from nexoclom.math import interpu
from pytest import approx
from nexoclom import __file__ as basefile

APLANET = 1.5
SPECIES = 'Na'

def test_gValue():
    """Compare gvalues with gvalues in IDL code."""
    ## Test 1
    g = gValue(SPECIES, 5891, APLANET)

    datafile = os.path.join(os.path.dirname(basefile), 'data', 'g-values', 
                                            'Na.D2.KillenXXXX.dat')
    aplanet_line = open(datafile).readline()
    aplanet = float(aplanet_line.split('=')[1].strip())
    fromfile = pd.read_csv(datafile, sep=':', skiprows=1)
    test_vel = fromfile.iloc[:,0].values
    test_g = fromfile.iloc[:,1].values * aplanet**2/APLANET**2
    s = np.argsort(test_vel)
    test_vel, test_g = test_vel[s], test_g[s]
    newg = np.interp(g.velocity.value, test_vel, test_g)

    assert g.wavelength == 5891.*u.AA, 'Wavelength failure'
    assert g.species == SPECIES, 'Species Failure'
    assert g.aplanet == APLANET*u.au, 'aplanet failure'
    assert g.g.value == approx(newg), 'g-values failure'

def test_radpresconst():
    rp_const = RadPresConst(SPECIES, APLANET)

    gfiles = ['Na.3303.Killen2009.dat', 'Na.D2.KillenXXXX.dat', 
              'Na.D1.KillenXXXX.dat']
    waves = [3303., 5891., 5897.]
    rp_const_new = np.zeros_like(rp_const.accel)
    for wave, gfile in zip(waves, gfiles):
        datafile = os.path.join(os.path.dirname(basefile), 'data', 'g-values', 
                                                gfile)
        aplanet_line = open(datafile).readline()
        aplanet = float(aplanet_line.split('=')[1].strip())

        fromfile = pd.read_csv(datafile, sep=':', skiprows=1)
        test_vel = fromfile.iloc[:,0].values * u.km/u.s
        test_g = fromfile.iloc[:,1].values * aplanet**2/APLANET**2 / u.s
        s = np.argsort(test_vel)
        test_vel, test_g = test_vel[s], test_g[s]

        # newg = np.interp(rp_const.velocity.value, test_vel, test_g)/u.s
        newg = interpu(rp_const.velocity, test_vel, test_g)
        rp_ = const.h/atomicmass(SPECIES)/(wave*u.AA) * newg
        rp = rp_.to(u.km/u.s**2)
        rp_const_new += rp

        # from IPython import embed; embed(); 

    assert rp_const.species == SPECIES
    assert rp_const.aplanet == APLANET*u.au
    assert rp_const.accel.value == approx(rp_const_new.value)
    assert rp_const.accel.unit is rp_const_new.unit

if __name__ == '__main__':
    test_gValue()
    test_radpresconst()
    
# import pickle
# g2, q_, q, rr = pickle.load(open('5891.0 Angstrom.pkl', 'rb'))
# print(np.all(np.isclose(g2, newg)))
# print(np.all(np.isclose(q_, rp_)))
# print(np.all(np.isclose(q, rp)))
# print(np.all(np.isclose(rr, rp_const_new)))
