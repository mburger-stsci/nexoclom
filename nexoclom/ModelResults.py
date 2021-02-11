import os.path
import numpy as np
import astropy.units as u
import mathMB
from atomicdataMB import gValue
from .database_connect import database_connect
from .input_classes import InputError
from .Output import Output


class ModelResult:
    def __init__(self, format, species=None):
        self.inputs = None
        if isinstance(format, str):
            if os.path.exists(format):
                self.format = {}
                with open(format, 'r') as f:
                    for line in f:
                        if ';' in line:
                            line = line[:line.find(';')]
                        elif '#' in line:
                            line = line[:line.find('#')]
                        else:
                            pass
                        
                        if '=' in line:
                            p, v = line.split('=')
                            self.format[p.strip().lower()] = v.strip()
                        else:
                            pass
            else:
                raise FileNotFoundError('ModelResult.__init__',
                                        'Format file not found.')
        elif isinstance(format, dict):
            self.format = format
        else:
            raise TypeError('ModelResult.__init__',
                            'format must be a dict or filename.')
            
        self.totalsource = 0.
        # Do some validation
        quantities = ['column', 'radiance', 'density']

        self.quantity = self.format.get('quantity', None)
        if (self.quantity is None) or (self.quantity not in quantities):
            raise InputError('ModelImage.__init__',
                             "quantity must be 'column' or 'radiance'")
        else:
            pass

        if self.quantity == 'radiance':
            # Note - only resonant scattering currently possible
            self.mechanism = ['resonant scattering']
    
            if 'wavelength' in self.format:
                self.wavelength = tuple(
                    int(m.strip())*u.AA
                    for m
                    in self.format['wavelength'].split(','))
            elif species is None:
                raise InputError('ModelImage.__init__',
                                 'Must provide either species or format.wavelength')
            elif species == 'Na':
                self.wavelength = (5891*u.AA, 5897*u.AA)
            elif species == 'Ca':
                self.wavelength = (4227*u.AA,)
            elif species == 'Mg':
                self.wavelength = (2852*u.AA,)
            else:
                raise InputError('ModelResult.__init__', ('Default wavelengths '
                              f'not available for {species}'))
        else:
            pass
    
    def search_for_outputs(self):
        self.outputfiles, self.outputpackets, self.outputsource = self.inputs.search()
        self.unit = u.def_unit('R_' + self.inputs.geometry.planet.object,
                               self.inputs.geometry.planet.radius)

    def transform_reference_frame(self, output):
        """If the image center is not the planet, transform to a
           moon-centric reference frame."""
        assert 0, 'Not ready yet.'

        # Load output

        # # Transform to moon-centric frame if necessary
        # if result.origin != result.inputs.geometry.planet:
        #     assert 0, 'Need to do transformation for a moon.'
        # else:
        #     origin = np.array([0., 0., 0.])*output.x.unit
        #     sc = 1.

        # Choose which packets to use
        # touse = output.frac >= 0 if keepall else output.frac > 0

        # packet positions relative to origin -- not rotated
        # pts_sun = np.array((output.x[touse]-origin[0],
        #                     output.y[touse]-origin[1],
        #                     output.z[touse]-origin[2]))*output.x.unit
        #
        # # Velocities relative to sun
        # vels_sun = np.array((output.vx[touse],
        #                      output.vy[touse],
        #                      output.vz[touse]))*output.vx.unit

        # Fractional content
        # frac = output.frac[touse]

        return output #, pts_sun, vels_sun, frac

    def packet_weighting(self, packets, out_of_shadow, aplanet):
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
                
            packets['weight'] = weight_resscat # + other stuff

        assert np.all(np.isfinite(packets['weight'])), 'Non-finite weights'
