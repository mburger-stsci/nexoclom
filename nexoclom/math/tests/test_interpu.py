import numpy as np
import astropy.units as u
from nexoclom.math import interpu
import pytest


def fn(x):
    return x**3 - 3*x**2 + 5

@pytest.mark.math
def test_interpu():
    xp = np.arange(20)
    fp = fn(xp)
    
    points = np.random.rand(30)*20-5
    result = np.interp(points, xp, fp)
    test = interpu(points*u.km, xp*u.km, fp*u.km**3)

    assert test.value == pytest.approx(result)
    assert test.unit == u.km**3
    
    test = interpu((points*u.km).to(u.imperial.mi), xp*u.km, fp*u.km**3)
    assert test.value == pytest.approx(result)

    with pytest.raises(TypeError):
        test = interpu(points, xp*u.km, fp*u.km**3)

    with pytest.raises(TypeError):
        _ = interpu(points*u.km, xp, fp*u.km**3)

    with pytest.raises(TypeError):
        _ = interpu(points*u.km, xp*u.km, fp)

    with pytest.raises(u.UnitConversionError):
        _ = interpu(points*u.s, xp*u.km, fp*u.km**3)
