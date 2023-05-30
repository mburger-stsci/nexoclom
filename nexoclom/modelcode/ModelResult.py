import os.path
import numpy as np
import copy
import astropy.units as u
from nexoclom import math as mathMB
from nexoclom.atomicdata import gValue
from nexoclom.modelcode.input_classes import InputError


class ModelResult:
    """Base class for nexoclom model comparisons with data.

    The ModelResult object is the base class for ModelImage (radiance and column
    density images), LOSResult (radiance or column along lines of sight), and
    ModelDensity (density along a trajectory or plane - not written yet).
    
    **Parameters**
    
    inputs
        An Input object
        
    params
        A dictionary containing parameters needed to create the model result.
        See LOSResult.py or ModelImage.py for details.
    
    **Methods**
    
    packet_weighting
        Determine the weighting factor each packet. When determining density
        or column, no weighting is needed. For radiance, takes into account
        determines g-value and takes into account shadowing.
        
    transform_reference_frame
        Will be used to transform packets from the central object reference frame
        to another location
        
    **Attributes**
    
    inputs
        Input object with the parameters associated with the model result
        
    params
        A dictionary containing parameters of the result. See LOSResult.py
        or ModelImage.py for details.
        
    outputfiles
        locations of saved Outputs associated with the inputs
    
    npackets
        total number of packets simulated
    
    totalsource
        total source in packets (if the initial fraction are all 1, then
        totalsource = npackets * nsteps
        
    quantity
        column, density, or radiance determined from params
        
    mechanism
        Emission mechanism if quantity = radiance else None
        
    wavelength
        Emssion wavelengths if quantity = radiance else None
      
    """
    def __init__(self, inputs, params):
        """
        :param inputs: Input object
        :param params: Dictionary with ModelResult parameters
        """
        self.inputs = copy.deepcopy(inputs)
        self.outid, self.outputfiles, _, _ = self.inputs.search()
        self.npackets = 0
        self.totalsource = 0.
        self.atoms_per_packet = 0.
        self.sourcerate = 0. * u.def_unit('10**23 atoms/s', 1e23 / u.s)
        if isinstance(params, str):
            if os.path.exists(params):
                self.params = {}
                with open(params, 'r') as f:
                    for line in f:
                        if ';' in line:
                            line = line[:line.find(';')]
                        elif '#' in line:
                            line = line[:line.find('#')]
                        else:
                            pass
                        
                        if '=' in line:
                            p, v = line.split('=')
                            self.params[p.strip().lower()] = v.strip()
                        else:
                            pass
            else:
                raise FileNotFoundError('ModelResult.__init__',
                                        'params file not found.')
        elif isinstance(params, dict):
            self.params = params
        else:
            raise TypeError('ModelResult.__init__',
                            'params must be a dict or filename.')
            
        # Do some validation
        quantities = ['column', 'radiance', 'density']
        self.quantity = self.params.get('quantity', None)
        if (self.quantity is None) or (self.quantity not in quantities):
            raise InputError('ModelImage.__init__',
                             "quantity must be 'column' or 'radiance'")
        else:
            pass

        if self.quantity == 'radiance':
            # Note - only resonant scattering currently possible
            self.mechanism = ['resonant scattering']
    
            if 'wavelength' in self.params:
                self.wavelength = tuple(sorted(int(m.strip())*u.AA for m
                    in self.params['wavelength'].split(',')))
            elif self.inputs.options.species is None:
                raise InputError('ModelImage.__init__',
                                 'Must provide either species or params.wavelength')
            elif self.inputs.options.species == 'Na':
                self.wavelength = (5891*u.AA, 5897*u.AA)
            elif self.inputs.options.species == 'Ca':
                self.wavelength = (4227*u.AA,)
            elif self.inputs.options.species == 'Mg':
                self.wavelength = (2852*u.AA,)
            else:
                raise InputError('ModelResult.__init__', ('Default wavelengths '
                                 f'not available for {self.inputs.options.species}'))
        else:
            self.mechanism = None
            self.wavelength = None
            
        self.unit = u.def_unit('R_' + self.inputs.geometry.planet.object,
                               self.inputs.geometry.planet.radius)
        
    def packet_weighting(self, packets, aplanet, out_of_shadow=1.):
        """
        Determine weighting factor for each packet
        :param packets: DataFrame with packet parameters
        :param out_of_shadow: Boolean array, True if in sunlight; False if in shadow
        :param aplanet: Distance of planet from Sun (used for g-value calculation)
        :return: Adds a 'weight' column to the packets DataFrame
        """
        if self.quantity == 'column':
            packets['weight'] = packets['frac']
        elif self.quantity == 'density':
            packets['weight'] = packets['frac']
        elif self.quantity == 'radiance':
            if 'resonant scattering' in self.mechanism:
                gg = np.zeros(len(packets))/u.s
                for w in self.wavelength:
                    gval = gValue(self.inputs.options.species, w, aplanet)
                    gg += mathMB.interpu(packets['radvel_sun'].values *
                                         self.unit/u.s, gval.velocity, gval.g)

                weight_resscat = packets['frac']*out_of_shadow*gg.value/1e6
            else:
                weight_resscat = np.zeros(len(packets))
                
            packets['weight'] = weight_resscat  # + other stuff
        else:
            raise InputError('ModelResults.packet_weighting',
                             f'{self.quantity} is invalid.')

        assert np.all(np.isfinite(packets['weight'])), 'Non-finite weights'
