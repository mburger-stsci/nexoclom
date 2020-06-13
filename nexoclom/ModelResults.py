import os.path
import numpy as np
import pandas as pd
import astropy.units as u
import mathMB
from atomicdataMB import gValue
from .database_connect import database_connect
from .input_classes import InputError
from .Output import Output


class ModelResult:
    def __init__(self, inputs, format, filenames=None, output=None):
        self.inputs = inputs
        if isinstance(output, Output):
            self.npackets = output.npackets
            self.totalsource = output.totalsource
            self.filenames = [None]
        elif isinstance(filenames, str):
            self.filenames = [filenames]
            with database_connect() as con:
                packs = pd.read_sql(f'''SELECT npackets, totalsource
                                        FROM outputfile
                                        WHERE filename={filenames}''', con)
            self.npackets = packs.npackets[0]
            self.totalsource = packs.totalsource[0]
        elif filenames is None:
            self.filenames, self.npackets, self.totalsource = inputs.search()
        else:
            raise Exception
            
        if self.npackets == 0:
            pass
        else:
            self.mod_rate = self.totalsource/inputs.options.endtime.value
            self.atoms_per_packet = 1e23/self.mod_rate
            print(f'Total number of packets run = {self.npackets}')
            print(f'Total source = {self.totalsource} packets')
            print(f'1 packet represents {self.atoms_per_packet} atoms')
            print(f'Model rate = {self.mod_rate} packets/sec')

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
            
            # Do some validation
            quantities = ['column', 'radiance', 'density']

            if 'quantity' in self.format:
                if self.format['quantity'] in quantities:
                    self.quantity = self.format['quantity']
                else:
                    raise InputError('ModelImage.__init__',
                                     "quantity must be 'column' or 'radiance'")
            else:
                raise InputError('ModelImage.__init__',
                                 'quantity must be specified.')

            if self.quantity == 'radiance':
                # Note - only resonant scattering currently possible
                self.mechanism = ['resonant scattering']
        
                if 'wavelength' in self.format:
                    self.wavelength = tuple(
                        int(m.strip())*u.AA
                        for m
                        in self.format['wavelength'].split(','))
                elif inputs.options.species == 'Na':
                    self.wavelength = (5891*u.AA, 5897*u.AA)
                elif inputs.options.species == 'Ca':
                    self.wavelength = (4227*u.AA,)
                elif inputs.options.species == 'Mg':
                    self.wavelength = (2852*u.AA,)
                else:
                    raise InputError('ModelResult.__init__', ('Default wavelengths '
                                  f'not available for {inputs.options.species}'))
            else:
                pass

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
