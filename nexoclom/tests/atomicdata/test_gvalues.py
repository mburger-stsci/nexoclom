import numpy as np
import pandas as pd
import astropy.units as u
from astropy import constants as const
from nexoclom.atomicdata import atomicmass, gValue, RadPresConst
from nexoclom.math import interpu
from pytest import approx
import pytest 

args_gval = [('Na', 5891, 1.5),
             ('Ca', 4227*u.AA, 0.3),
             ('X', 3333*u.AA, 1.0*u.au)]
args_rp = [('Na', 1.5),
           ('Ca', 0.3*u.au),
           ('X', 1.0*u.au)]

@pytest.mark.atomicdata
@pytest.mark.parametrize('species, wavelength, aplanet', args_gval)
def test_gValue(species, wavelength, aplanet):
    """Compare gvalues with gvalues in IDL code."""
    ## Test 1
    g = gValue(species, wavelength, aplanet)
    if g.filename is not None:
        aplan_line = open(g.filename).readline()
        aplan = float(aplan_line.split('=')[1].strip())
        fromfile = pd.read_csv(g.filename, sep=':', skiprows=1)
        test_vel = fromfile.iloc[:,0].values
        test_g = fromfile.iloc[:,1].values * aplan**2/aplanet**2
        s = np.argsort(test_vel)
        test_vel, test_g = test_vel[s], test_g[s]
        newg = np.interp(g.velocity.value, test_vel, test_g)

        assert g.species == species, 'Species Failure'

        if isinstance(wavelength, type(1*u.AA)):
            assert g.wavelength == wavelength, 'Wavelength failure'
        else:
            assert g.wavelength == wavelength*u.AA

        if isinstance(aplanet, type(1*u.au)):
            assert g.aplanet == aplanet, 'aplanet failure'
        else:
            assert g.aplanet == aplanet*u.au, 'aplanet failure'

        assert g.g.value == approx(newg), 'g-values failure'

@pytest.mark.atomicdata
@pytest.mark.parametrize('species, aplanet', args_rp)
def test_radpresconst(species, aplanet):
    rp_const = RadPresConst(species, aplanet)

    rp_const_new = np.zeros_like(rp_const.accel)
    for wave, in rp_const.wavelength:
        g = gValue(species, wave, aplanet)
        aplan_line = open(g.filename).readline()
        aplan = float(aplan_line.split('=')[1].strip())

        fromfile = pd.read_csv(g.filename, sep=':', skiprows=1)
        test_vel = fromfile.iloc[:,0].values * u.km/u.s
        test_g = fromfile.iloc[:,1].values * aplan**2/aplanet**2 / u.s
        s = np.argsort(test_vel)
        test_vel, test_g = test_vel[s], test_g[s]

        newg = interpu(rp_const.velocity, test_vel, test_g)
        rp_ = const.h/atomicmass(species)/(wave*u.AA) * newg
        rp = rp_.to(u.km/u.s**2)
        rp_const_new += rp

    assert rp_const.species == species

    if isinstance(aplanet, type(1*u.au)):
        assert rp_const.aplanet == aplanet
    else:
        assert rp_const.aplanet == aplanet*u.au
    assert rp_const.accel.value == approx(rp_const_new.value)
    assert rp_const.accel.unit is rp_const_new.unit

if __name__ == '__main__':
    test_gValue()
    test_radpresconst()