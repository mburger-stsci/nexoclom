import numpy as np
import astropy.units as u
from nexoclom.solarsystem import planet_dist, SSObject
import pytest

@pytest.mark.solarsystem
def test_planet_dist_circle():
    planet = SSObject('Mercury')
    planet.e = 0
    taa = np.arange(0, 2*np.pi, 100)*u.rad
    distance, radvel = planet_dist(planet, taa=taa)
    assert np.all(distance == planet.a), 'Radius should be constant for eps=0'
    assert np.all(radvel == 0*u.km/u.s), 'Radial velocity should be 0 for eps=0'
    
@pytest.mark.solarsystem
def test_plot_inputs():
    assert planet_dist(SSObject('Nothing')) is None
    assert planet_dist(SSObject('Mercury'), taa=None, time=None) is None
