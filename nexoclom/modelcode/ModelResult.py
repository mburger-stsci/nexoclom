import os.path
import sys
import numpy as np
import copy
import shutil
import astropy.units as u
import pandas as pd
import pickle
from sklearn.neighbors import BallTree
import time
from nexoclom import math as mathMB
from nexoclom.atomicdata import gValue
from nexoclom.modelcode.input_classes import InputError
from nexoclom.modelcode.Output import Output
from nexoclom.modelcode.SourceMap import SourceMap
from nexoclom import __file__ as basefile
try:
    import condorMB
except:
    pass


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

    def make_source_map(self, smear_radius=10*np.pi/180, nlonbins=180,
                        nlatbins=90, nvelbins=100, nazbins=90, naltbins=23,
                        use_condor=False, normalize=True):
        """
        At each point in lon/lat grid want:
            * Source flux (atoms/cm2/s
            * Speed distribution (f_v vs v)
            * Azimuthal distribution (f_az vs az) -> measured CCW from north
            * Altitude distribution (f_alt vs alt) -> tangent = 0, normal = 90
        """
        source, available = {}, {}
        for outputfile in self.outputfiles:
            output = Output.restore(outputfile)
            included = output.X.loc[output.X.frac > 0, 'Index'].unique()
            X0 = output.X0[['x', 'y', 'z', 'vx', 'vy', 'vz', 'frac']]
            X0.loc[included, 'included'] = True
            del output
            
            X0.included.fillna(False, inplace=True)
            velocity = (X0[['vx', 'vy', 'vz']].values *
                        self.inputs.geometry.planet.radius.value)
            speed = np.linalg.norm(velocity, axis=1)

            # Radial, east, north unit vectors
            rad = X0[['x', 'y', 'z']].values
            rad_ = np.linalg.norm(rad, axis=1)
            rad = rad/rad_[:, np.newaxis]

            east = X0[['y', 'x', 'z']].values
            east[:,1] = -1*east[:,1]
            east[:,2] = 0
            east_ = np.linalg.norm(east, axis=1)
            east = east/east_[:, np.newaxis]

            # north = np.array([-X0.z.values * X0.x.values,
            #                   -X0.z.values * X0.y.values,
            #                   X0.x.values**2 + X0.y.values**2])
            north = np.cross(rad, east, axis=1)
            north_ = np.linalg.norm(north, axis=1)
            north = north/north_[:, np.newaxis]

            v_rad = np.sum(velocity * rad, axis=1)
            v_east = np.sum(velocity * east, axis=1)
            v_north = np.sum(velocity * north, axis=1)

            v_rad_over_speed = v_rad/speed
            v_rad_over_speed[v_rad_over_speed > 1] = 1
            v_rad_over_speed[v_rad_over_speed < -1] = -1

            assert np.all(np.isclose(v_rad**2 + v_east**2 + v_north**2, speed**2))
            X0.loc[:, 'altitude'] = np.arcsin(v_rad_over_speed)
            X0.loc[:, 'azimuth'] = (np.arctan2(v_north, v_east) + 2*np.pi) % (2*np.pi)
            X0.loc[:, 'v_rad'] = v_rad
            X0.loc[:, 'v_east'] = v_east
            X0.loc[:, 'v_north'] = v_north
            X0.loc[:, 'speed'] = speed
            X0.loc[:, 'longitude'] = (np.arctan2(X0.x.values, -X0.y.values) + 2*np.pi) % (2*np.pi)
            X0.loc[:, 'latitude'] = np.arcsin(X0.z.values)
            
            
            tree = BallTree(X0[['latitude', 'longitude']], metric='haversine')

            X0_ = X0[X0.frac > 0]
            # Calculate the histograms and available fraction
            for which in (0, 1):
                if which == 0:
                    distribution = copy.deepcopy(source)
                    weight = X0_.frac
                else:
                    distribution = copy.deepcopy(available)
                    weight = np.ones_like(X0_.frac.values)

                if (nlonbins > 0) and (nlatbins > 0):
                    abundance = mathMB.Histogram2d(X0_.longitude, X0_.latitude, weights=weight,
                                                range=[[0, 2*np.pi], [-np.pi/2, np.pi/2]],
                                                bins=(nlonbins, nlatbins))
                    abundance.x, abundance.dx = abundance.x * u.rad, abundance.dx * u.rad
                    abundance.y, abundance.dy = abundance.y * u.rad, abundance.dy * u.rad

                    gridlatitude, gridlongitude = np.meshgrid(abundance.y, abundance.x)

                    points = np.array([gridlatitude.flatten(), gridlongitude.flatten()]).T
                    ind = tree.query_radius(points, smear_radius)
                    fraction_observed = np.ndarray((points.shape[0], ))
                    for index in range(points.shape[0]):
                        included = X0.loc[ind[index], 'included']
                        fraction_observed[index] = sum(included)/included.shape[0]

                    if 'abundance' in distribution:
                        distribution['abundance'] += abundance.histogram
                        distribution['fraction_observed'] += fraction_observed.reshape(
                            gridlongitude.shape)/len(self.outputfiles)
                    else:
                        distribution['abundance'] = abundance.histogram
                        distribution['longitude'] = abundance.x
                        distribution['latitude'] = abundance.y
                        distribution['fraction_observed'] = fraction_observed.reshape(
                            gridlongitude.shape)/len(self.outputfiles)
                else:
                    pass

                if nvelbins > 0:
                    velocity = mathMB.Histogram(X0_.speed, bins=nvelbins,
                                                range=[0, X0_.speed.max()],
                                                weights=weight)
                    velocity.x = velocity.x * u.km / u.s
                    velocity.dx = velocity.dx * u.km / u.s
                    
                    if 'speed' in distribution:
                        distribution['speed_dist'] += velocity.histogram
                    else:
                        distribution['speed'] = velocity.x
                        distribution['speed_dist'] = velocity.histogram
                else:
                    pass

                if naltbins > 0:
                    altitude = mathMB.Histogram(X0_.altitude, bins=naltbins,
                                                range=[0, np.pi / 2], weights=weight)
                    altitude.x = altitude.x * u.rad
                    altitude.dx = altitude.dx * u.rad
                    
                    if 'altitude' in distribution:
                        distribution['altitude_dist'] += altitude.histogram
                    else:
                        distribution['altitude'] = altitude.x
                        distribution['altitude_dist'] = altitude.histogram
                else:
                     pass
                
                if nazbins > 0:
                    azimuth = mathMB.Histogram(X0_.azimuth, bins=nazbins,
                                               range=[0, 2 * np.pi], weights=weight)
                    azimuth.x = azimuth.x * u.rad
                    azimuth.dx = azimuth.dx * u.rad
                    
                    if 'azimuth' in distribution:
                        distribution['azimuth_dist'] += azimuth.histogram
                    else:
                        distribution['azimuth'] = azimuth.x
                        distribution['azimuth_dist'] = azimuth.histogram
                else:
                    pass
                    
            if which == 0:
                source = distribution
            else:
                available = distribution
                
            del X0, X0_, output
            
        ## normalization
        if normalize:
            for which in (0, 1):
                if which == 0:
                    distribution = copy.deepcopy(source)
                else:
                    distribution = copy.deepcopy(available)
                    
                if 'abundance' in distribution:
                    # Convert histogram to flux
                    # (a) divide by area of a grid cell
                    #   Surface area of a grid cell =
                    #       R**2 (lambda_2 - lambda_1) (sin(phi2)-sin(phi1))
                    #   https://onlinelibrary.wiley.com/doi/epdf/10.1111/tgis.12636, eqn 1
                    # (b) Multiply by source rate
                    dx = distribution['longitude'][1] - distribution['longitude'][0]
                    dy = distribution['latitude'][1] - distribution['latitude'][0]
                    _, gridlatitude = np.meshgrid(distribution['longitude'],
                                                  distribution['latitude'])
                    area = (self.inputs.geometry.planet.radius**2 * dx.value *
                            (np.sin(gridlatitude + dy / 2) -
                             np.sin(gridlatitude - dy / 2)))
        
                    source.histogram = (distribution['abundance'] /
                                        distribution['abundance'].sum() /
                                        area.T.to(u.cm**2) * self.sourcerate.to(1 / u.s))
                    from IPython import embed; embed()
                    import sys; sys.exit()
                else:
                    pass
                
            if 'speed' in distribution:
                dv = distribution['speeed'][1] - distribution['speed'][0]
                distribution['speed_dist'] = (self.sourcerate *
                                              distribution['speed_dist'] /
                                              distribution['speed_dist'].sum() / dv)
                distribution['speed_dist'] = distribution['speed_dist'] * (
                    self.sourcerate.unit * u.def_unit('(km/s)^-1', u.s/u.km))
            else:
                pass
            
            if 'altitude' in distribution:
                dalt = distribution['altitude'][1] - distribution['altitude'][0]
                distribution['altitude'] = (self.sourcerate * distribution['altitude'] /
                                            distribution['altitude'].sum() / dalt)
            else:
                pass

            if 'azimuth' in distribution:
                daz = distribution['azimuth'][1] - distribution['azimuth'][0]
                distribution['azimuth'] = (self.sourcerate * distribution['azimuth'] /
                                            distribution['azimuth'].sum() / daz)
            else:
                pass
            
            if which == 0:
                source = distribution
            else:
                available = distribution
        else:
            pass
            
        return source, available
