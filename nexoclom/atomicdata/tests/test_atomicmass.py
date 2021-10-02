import astropy.units as u
import periodictable as pt
import random
from nexoclom.atomicdata import atomicmass
import pytest


sp = [random.choice(list(pt.elements)).symbol for _ in range(3)]
sp.extend(['H2O', 'CO2', 'C2H2OH'])
res = list([pt.formula(s).mass for s in sp])
sp.append('X')
res.append(None)

args = [pytest.param(s, r, id=s) for s,r in zip(sp, res)]

@pytest.mark.atomicdata
@pytest.mark.parametrize('species, result', args)
def test_atomicmass(species, result):
    if result is None:
        assert atomicmass(species) is None
    else:
        assert atomicmass(species) == result * u.u, f'{species} mass failure'