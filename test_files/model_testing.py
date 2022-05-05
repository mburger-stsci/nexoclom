import numpy as np
import os
import astropy.units as u
from nexoclom import Input
from MESSENGERuvvs import MESSENGERdata


orbit = 3576
species = 'Ca'
overwrite = True
npack = 5e5

# Load model results
mdata = MESSENGERdata(species, f'orbit={orbit}', load_spectra=False)
# mdata = MESSENGERdata(species, f'TAA >= 0 and TAA < 0.03489', load_spectra=False)
inputpath = '/Users/mburger/Work/Research/Mercury/inputfiles'

inputfile = os.path.join(inputpath, f'{species}.isotropic.flat.input')
inputs = Input(inputfile)
inputs.geometry.taa = np.median(mdata.data.taa)*u.rad
inputs.options.fitted = False

minalt = {'Ca':100, 'Na':50}
masking = f'minalt{minalt[species]};minsnr2'
# masking = f'minalt{[species]};minsnr2;siglimit10'
# masking = f'minsnr2'

mdata.model(inputs, npack, masking=masking, dphi=1 * u.deg,
            label='Uniform Model', overwrite=overwrite, packs_per_it=100000)

inputs.options.fitted = True
mdata.model(inputs, npack, masking=masking, label='Fitted Model, packs_per_it=100000')

mdata.plot(plot_method='bokeh', show=True,
           filename=f'{species}_Oribt{orbit:04d}.html')
