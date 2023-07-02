""" Tests whether a simple model will run to completion
"""
import os
import pytest
from nexoclom import Input, __path__
from MESSENGERuvvs import MESSENGERdata


inputfiles = ('Ca.reference.input', 'Na.reference.input')
@pytest.mark.integration
@pytest.mark.parametrize('inputfile', inputfiles)
def test_run_through(inputfile):
    basepath = os.path.dirname(__path__[0])
    inputfile_ = os.path.join(basepath, 'tests', 'test_data', 'inputfiles',
                              inputfile)
    inputs = Input(inputfile_)
    species = inputfile[:2]
    print(species)
    mdata = MESSENGERdata(species, 'orbit=36')

    npack, pperit = 2e4, 1e4
    minalt = {'Ca':100, 'Na':50, 'Mg':50}
    mask = f'minalt{minalt[species]};minsnr2'
    
    # Test with delay
    mdata.model(inputs, npack, masking=mask, label='Test', overwrite=True,
                packs_per_it=pperit, distribute='delay')
    
    # Test with serial
    mdata.model(inputs, npack, masking=mask, label='Test', overwrite=True,
                packs_per_it=pperit, distribute=None)


if __name__ == '__main__':
    for inputfile in inputfiles:
        test_run_through(inputfile)
