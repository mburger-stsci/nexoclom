import numpy as np
from ..LossInfo import LossInfo


def test_LossInfo():
    rate0 = LossInfo('Na', -2., 0.3)
    assert rate0.photo == 0.5
    
    rate1 = LossInfo('Ca', 0, 1)
    assert rate1.photo == 7e-5
    
    rate2 = LossInfo('Mg', 0, 0.3)
    assert np.isclose(rate2.photo, 7.21111111111e-6)