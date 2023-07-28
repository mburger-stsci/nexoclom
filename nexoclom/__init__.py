from importlib.metadata import version
from nexoclom.utilities.configure import configure
config, engine = configure(verbose=False)

import dask
# dask.config.set(scheduler='processes')
# from dask.distributed import Client

from nexoclom.initial_state.Input import Input
from nexoclom.particle_tracking.Output import Output
from nexoclom.data_simulation.LOSResult import LOSResult
from nexoclom.data_simulation.LOSResultFitted import LOSResultFitted
from nexoclom.data_simulation.ModelImage import ModelImage
from nexoclom.solarsystem import SSObject


__name__ = 'nexoclom'
__author__ = 'Matthew Burger'
__email__ = 'mburger@stsci.edu'
__version__ = version("nexoclom")
__date__ = '2023-07-28'
