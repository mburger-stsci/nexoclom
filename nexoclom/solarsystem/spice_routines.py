""" Wrappers for working with spice
SpiceyPy reference: Annex et al., (2020). SpiceyPy: a Pythonic Wrapper for the
SPICE Toolkit. Journal of Open Source Software, 5(46), 2050,
https://doi.org/10.21105/joss.02050
"""
import os
import spiceypy as spice
import requests
from nexoclom import config
from nexoclom.utilities import ConfigfileError


def load_kernels():
    """Load generic spice kernels. Retrieves if necessary."""
    # leap second kernel
    datapath = config.__dict__.get('datapath', None)
    if datapath is None:
        raise ConfigfileError('datapath')
    else:
        pass
        
    baseurl = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/'
    kernel_path = os.path.join(datapath, 'SpiceKernels')
    lsk_kernel = os.path.join(kernel_path, 'lsk', 'naif0012.tls')
    if not os.path.exists(lsk_kernel):
        if not os.path.exists(os.path.dirname(lsk_kernel)):
            os.makedirs(os.path.dirname(lsk_kernel))
        lsk_url = f'{baseurl}/lsk/latest_leapseconds.tls'
        print(f'Retreiving leapsecond kernel {os.path.basename(lsk_kernel)}')
        page = requests.get(lsk_url).text
        with open(lsk_kernel, 'w') as file:
            file.write(page)
    else:
        pass
    
    gm_kernel = os.path.join(kernel_path, 'pck', 'gm_de440.tpc')
    if not os.path.exists(gm_kernel):
        if not os.path.exists(os.path.dirname(gm_kernel)):
            os.makedirs(os.path.dirname(gm_kernel))
        gm_url = f'{baseurl}/pck/{os.path.basename(gm_kernel)}'
        print(f'Retreiving GM kernel {os.path.basename(gm_kernel)}')
        page = requests.get(gm_url).text
        with open(gm_kernel, 'w') as file:
            file.write(page)
    else:
        pass

    pck_kernel = os.path.join(kernel_path, 'pck', 'pck00011.tpc')
    if not os.path.exists(pck_kernel):
        if not os.path.exists(os.path.dirname(pck_kernel)):
            os.makedirs(os.path.dirname(pck_kernel))
        pck_url = f'{baseurl}/pck/{os.path.basename(pck_kernel)}'
        print(f'Retreiving pck kernel {os.path.basename(pck_kernel)}')
        page = requests.get(pck_url).text
        with open(pck_kernel, 'w') as file:
            file.write(page)
    else:
        pass

    spice.furnsh(lsk_kernel)
    spice.furnsh(gm_kernel)
    spice.furnsh(pck_kernel)

if __name__ == '__main__':
    load_kernels()
    from inspect import currentframe, getframeinfo
    frameinfo = getframeinfo(currentframe())
    print(frameinfo.filename, frameinfo.lineno)
    from IPython import embed; embed()
    import sys; sys.exit()
