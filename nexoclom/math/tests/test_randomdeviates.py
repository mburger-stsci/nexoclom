import numpy as np
from nexoclom.math import random_deviates_1d, random_deviates_2d
import nexoclom.math.distributions as distributions
import pytest



@pytest.mark.math
def test_random_deviates():
    # Test against gaussian
    a, b, c = list(np.random.random(3)*5)
    x0 = np.linspace(-10, 10, 1000)
    f0 = a * np.exp((x0 - b)**2/c**2)

    dev = random_deviates_1d(x0, f0, 1e5)
    
    from IPython import embed; embed()
    import sys; sys.exit()
    