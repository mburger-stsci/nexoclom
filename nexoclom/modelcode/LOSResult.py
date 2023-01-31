import numpy as np
import pandas as pd
from sklearn.neighbors import BallTree
import pickle
import astropy.units as u
from astropy.modeling import models, fitting
from astropy.visualization import PercentileInterval
import sqlalchemy as sqla
import dask

import nexoclom.math as mathMB
from nexoclom import engine
from nexoclom.modelcode.Output import Output
from nexoclom.modelcode.SourceMap import SourceMap
from nexoclom.modelcode.ModelResult import ModelResult
from nexoclom.modelcode.compute_iteration import (compute_iteration,
                                                  IterationResult,
                                                  IterationResultFitted)
from nexoclom import __file__ as basefile


class LOSResult(ModelResult):
    """Class to contain the LOS result from multiple outputfiles.
    
    Determine column or emission along lines of sight.
    This assumes the model has already been run.

    **Parameters**
    
    scdata
        Spacecraft data object (currently designed for MESSENGERdata object
        but can be faked for other types of data)

    params
        A dictionary containing the keys
        
            * quantity [required]: column, density, radiance
            
            * wavelength [optional]: For radiance, wavelenghts to be simulated.
            If not given, uses defaults for species. Must be a valid emission
            line for the species.
            
        More parameters will be added when more emission processes are included.
        For now, the easiest is `params = {'format': 'radiance'}`

    dphi
        Angular size of the view cone. Default = r deg.
        
    **Methods**
    
    **Attributes**
   
    species, query
        The species and query used to retrieve the data used. These can be
        used to retrieve the data if necessary
        
    type
        'LineOfSight' for a line of sight result
        
    dphi
        boresight opening angle
        
    radiance
        Pandas series containing modeled radiance along each line of sight
        
    npackets
        Pandas series containing the number of packets along each line of sight
    
    sourcemap
        Characterization of the initial source (spatial and velocity distributions)
        
    modelfiles
        Saved LOS Iteration results
    
    _oedge
        Maximum distance from the s/c to integrate. Twice the outer edge of the
        simulation region or 100 R_planet, whichever is less.
    """
    def __init__(self, scdata, inputs, params=None, dphi=1*u.deg, **kwargs):
        """Initializes the LOSResult and runs the model if necessary"""
        if params is None:
            params = {'quantity': 'radiance'}
        else:
            pass

        scdata.set_frame('Model')
        super().__init__(inputs, params)
        
        # Basic information
        self.species = scdata.species
        self.query = scdata.query
        self.type = 'LineOfSight'
        self.dphi = dphi.to(u.rad).value
        self._oedge = np.min([self.inputs.options.outeredge*2, 100])

        self.fitted = self.inputs.options.fitted
        nspec = len(scdata)
        self.radiance = pd.Series(np.zeros(nspec), index=scdata.data.index)
        self.radiance_unit = u.def_unit('kR', 1e3*u.R)
        self.sourcemap = None
        self.modelfiles = None
        
        self.goodness_of_fit = None
        self.mask = None
        self.masking = kwargs.get('masking', None)
        self.fit_method = kwargs.get('fit_method', None)
        self.label = kwargs.get('label', 'LOSResult')

    def __repr__(self):
        return self.__str__()
        
    def __str__(self):
        return f'''Model Label = {self.label}
quantity = {self.quantity}
npackets = {self.npackets}
totalsource = {self.totalsource}
atoms per packet = {self.atoms_per_packet}
sourcerate = {self.sourcerate}
dphi = {self.dphi}
fit_method = {self.fit_method}
fitted = {self.fitted}'''
    
    def search_iterations(self, fitted=False):
        """
        :return: dictionary containing search results:
                 {outputfilename: (modelfile_id, modelfile_name)}
        """
        search_results = {}
        for oid, outputfile in zip(self.outid, self.outputfiles):
            metadata_obj = sqla.MetaData()
            table = sqla.Table("uvvsmodels", metadata_obj, autoload_with=engine)

            ufit_id = (self.unfit_outid if fitted else None)
            query = sqla.select(table).where(
                table.columns.out_idnum == oid,
                table.columns.unfit_idnum == ufit_id,
                table.columns.quantity == self.quantity,
                table.columns.query == self.query,
                table.columns.dphi == self.dphi,
                table.columns.mechanism == self.mechanism,
                table.columns.wavelength == [w.value for w in self.wavelength],
                table.columns.fitted == fitted)
                
            with engine.connect() as con:
                result = pd.DataFrame(con.execute(query))
                
            # Should only have one match per outputfile
            assert len(result) <= 1
            
            if len(result) == 0:
                search_results[outputfile] = None
            else:
                search_results[outputfile] = (result.loc[0, 'idnum'],
                                              result.loc[0, 'unfit_idnum'],
                                              result.loc[0, 'filename'])
        
        return search_results
    
    def restore_iteration(self, search_result, save_ufit_id=False):
        # Restore is on an outputfile basis
        idnum, ufit_idnum, modelfile = search_result
        print(f'Restoring modelfile {modelfile}.')
        with open(modelfile, 'rb') as f:
            iteration_result = pickle.load(f)
        
        iteration_result.modelfile = modelfile
        iteration_result.model_idnum = idnum
        if save_ufit_id:
            self.ufit_idnum = ufit_idnum
        else:
            pass
        
        return iteration_result
    
    def make_mask(self, data):
        mask = np.array([True for _ in data.radiance])
        sigmalimit = None
        if self.masking is not None:
            for masktype in self.masking.split(';'):
                masktype = masktype.strip().lower()
                if masktype.startswith('middle'):
                    perinterval = float(masktype[6:])
                    # Estimate model strength (source rate) by fitting middle %
                    interval = PercentileInterval(perinterval)
                    lim = interval.get_limits(data)
                    mask = (mask &
                            (data.radiance >= lim[0]) &
                            (data.radiance <= lim[1]))
                elif masktype.startswith('minalt'):
                    minalt = float(masktype[6:])
                    mask = mask & (data.alttan >= minalt)
                elif masktype.startswith('minsnr'):
                    minSNR = float(masktype[6:])
                    snr = data.radiance/data.sigma
                    mask = mask & (snr > minSNR)
                elif masktype.startswith('siglimit'):
                    sigmalimit = float(masktype[8:])
                else:
                    raise ValueError('nexoclom.math.fit_model',
                                     f'masking = {masktype} not defined.')
        else:
            pass
        
        return mask, sigmalimit

    def simulate_data_from_inputs(self, scdata, distribute=None):
        """Given a set of inputs, determine what the spacecraft should see.
        Models should have already been run.
        
        **Outputs**
        """
        # If using a planet-fixed source map, need to set subsolarlon
        if ((self.inputs.spatialdist.type == 'surface map') and
            (self.inputs.spatialdist.coordinate_system == 'planet-fixed')):
            self.inputs.spatialdist.subsolarlon = scdata.subslong.median() * u.rad
        else:
            pass
    
        # Find the output files that have already been run for these inputs
        (self.outid, self.outputfiles, self.npackets,
         self.totalsource) = self.inputs.search()
        print(f'LOSResult: {len(self.outid)} output files found.')
        if self.npackets == 0:
            raise RuntimeError('No packets found for these Inputs.')
    
        # Find any model results that have been run for these inputs
        data = scdata.data
        search_results = self.search_iterations()
        
        while None in search_results.values():
            # Will retry if something fails due to memory error
            print(f'LOSResult: {list(search_results.values()).count(None)} '
                  'to compute')
            if distribute in ('delay', 'delayed'):
                iterations = [dask.delayed(compute_iteration)(self, outputfile,
                                                              scdata, True)
                              for outputfile, search_result in search_results.items()
                              if search_result is None]
                dask.compute(*iterations)
            else:
                for outputfile, search_result in search_results.items():
                    if search_result is None:
                        compute_iteration(self, outputfile, scdata)
                    else:
                        pass
                    
            search_results = self.search_iterations()
            
        iteration_results = []
        for outputfile, search_result in search_results.items():
            assert search_result is not None
            iteration_result = self.restore_iteration(search_result)
            iteration_result.model_idnum = search_result[0]
            iteration_result.modelfile = search_result[2]
            assert len(iteration_result.radiance) == len(data)
            iteration_results.append(iteration_result)
        else:
            pass
    
        # combine iteration_results
        self.modelfiles = {}
        for iteration_result in iteration_results:
            self.radiance += iteration_result.radiance
            self.modelfiles[iteration_result.outputfile] = iteration_result.modelfile
    
        # need model rate for this output
        model_rate = self.totalsource / self.inputs.options.endtime.value
        self.atoms_per_packet = 1e23 / model_rate
        self.radiance *= self.atoms_per_packet/1e3  # kR
        self.determine_source_rate(scdata)
        self.atoms_per_packet *= self.sourcerate.unit * u.s
        self.outputfiles = list(self.modelfiles.keys())
    
        print(self.totalsource, self.atoms_per_packet)
        
    def determine_source_rate(self, scdata):
        mask, sigmalimit = self.make_mask(scdata.data)
        linmodel = models.Multiply()
        fitter = fitting.LinearLSQFitter()
        best_fit = fitter(linmodel, self.radiance.values[mask],
                          scdata.data.radiance.values[mask],
                          weights=1./scdata.data.sigma.values[mask]**2)
        
        if sigmalimit is not None:
            diff = np.abs((scdata.data.radiance.values -
                           best_fit.factor*self.radiance.values)/
                          scdata.data.sigma)
            mask = mask & (diff < sigmalimit)
            best_fit = fitter(linmodel, self.radiance.values[mask],
                              scdata.data.radiance.values[mask],
                              weights=1./scdata.data.sigma.values[mask]**2)
        else:
            pass

        self.radiance *= best_fit.factor.value
        self.sourcerate = best_fit.factor.value * u.def_unit('10**23 atoms/s', 1e23 / u.s)
        self.goodness_of_fit = None

        self.mask = mask

    def make_source_map(self, smear_radius=np.radians(10), nlonbins=180,
                        nlatbins=90, nvelbins=100, nazbins=90, naltbins=23,
                        use_condor=False, normalize=True, do_source=True,
                        do_available=True):
        """
        At each point in lon/lat grid want:
            * Source flux (atoms/cm2/s
            * Speed distribution (f_v vs v)
            * Azimuthal distribution (f_az vs az) -> measured CCW from north
            * Altitude distribution (f_alt vs alt) -> tangent = 0, normal = 90
        """
        todo, source, available = [], None, None
        if do_source:
            todo.append(0)
            source = {}
        else:
            pass

        if do_available:
            todo.append(1)
            available = {}
        else:
            pass
        
        for which in todo:
            if which == 0:
                print('Determining modeled source')
            else:
                print('Determining available source')
            distribution = {'abundance_uncor': np.zeros((nlonbins, nlatbins)) / u.s,
                            'abundance': np.zeros((nlonbins, nlatbins)) / u.s,
                            'fraction_observed': np.zeros((nlonbins, nlatbins)),
                            'speed_dist_map': np.zeros((nlonbins, nlatbins,
                                                        nvelbins)) * u.s/u.km,
                            'n_included': np.zeros((nlonbins, nlatbins)),
                            'n_total': np.zeros((nlonbins, nlatbins)),
                            'speed_dist': np.zeros(nvelbins) * u.s / u.km,
                            'altitude_dist': np.zeros(naltbins) / u.rad,
                            'azimuth_dist': np.zeros(nazbins) / u.rad}
            vmax = None

            for outputfile, modelfile in self.modelfiles.items():
                print(outputfile)
                output = Output.restore(outputfile)
                
                if output.inputs.options.fitted:
                    with open(modelfile, 'rb') as file:
                        model = pickle.load(file)
                    
                    with open(model.unfit_modelfile, 'rb') as file:
                        unfit_model = pickle.load(file)
                    
                    included = set()
                    for row in unfit_model.used_packets0:
                        included = included.union(row)
                    included = list(included)
                else:
                    included = list(output.X0.index)
                
                X0 = output.X0
                X0.loc[included, 'included'] = True
                X0.included.fillna(False, inplace=True)
                del output
                
                velocity = (X0[['vx', 'vy', 'vz']].values *
                            self.inputs.geometry.planet.radius.value)
                speed = np.linalg.norm(velocity, axis=1)
                if vmax is None:
                    vmax = np.ceil(speed.max())
                else:
                    pass

                # Convert v_x, v_y, v_z -> v_east, v_north, v_radial vectors
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

                assert np.all(np.isclose(v_rad**2 + v_east**2 + v_north**2,
                                         speed**2))
                X0.loc[:, 'altitude'] = np.arcsin(v_rad_over_speed)
                X0.loc[:, 'azimuth'] = (np.arctan2(v_north, v_east) +
                                        2*np.pi) % (2*np.pi)
                X0.loc[:, 'v_rad'] = v_rad
                X0.loc[:, 'v_east'] = v_east
                X0.loc[:, 'v_north'] = v_north
                X0.loc[:, 'speed'] = speed
                # X0.loc[:, 'longitude'] = (np.arctan2(X0.x.values, -X0.y.values) + 2*np.pi) \
                #                          % (2*np.pi)
                # X0.loc[:, 'latitude'] = np.arcsin(X0.z.values)
                tree = BallTree(X0[['latitude', 'longitude']], metric='haversine')
                
                X0_subset = X0[X0.included]
                if which == 0:
                    weight = X0.frac
                    weight_subset = X0_subset.frac
                else:
                    weight = np.ones_like(X0.frac.values)
                    weight_subset = np.ones_like(X0_subset.frac.values)
                
                # Calculate the histograms and available fraction
                abundance = mathMB.Histogram2d(X0_subset.longitude,
                                               X0_subset.latitude,
                                               weights=weight_subset,
                                               range=[[0, 2*np.pi],
                                                      [-np.pi/2, np.pi/2]],
                                               bins=(nlonbins, nlatbins))
                gridlatitude, gridlongitude = np.meshgrid(abundance.y,
                                                          abundance.x)

                points = np.array([gridlatitude.flatten(),
                                   gridlongitude.flatten()]).T
                ind = tree.query_radius(points, smear_radius)
                n_included = np.zeros((points.shape[0], ))
                n_total = np.zeros((points.shape[0], ))
                v_point = np.zeros((points.shape[0], nvelbins))
                for index in range(points.shape[0]):
                    if len(ind[index]) > 0:
                        included = X0.loc[ind[index], 'included']
                        weight_ = X0.loc[ind[index], 'frac']
                        n_included[index] = np.sum(included*weight_)/np.sum(weight_)
                        n_total[index] = weight_.sum()
                        
                        # No weighting because assumption is all the atoms are ejected
                        # from the same point (points[index, :])
                        # n_included[index] = np.sum(included)
                        # n_total[index] = len(included)
                        vpoint_ = mathMB.Histogram(X0.loc[ind[index], 'speed'],
                                                   bins=nvelbins, range=[0, vmax])
                        v_point[index,:] += vpoint_.histogram
                    else:
                        pass

                distribution['abundance_uncor'] += abundance.histogram / u.s
                distribution['longitude'] = abundance.x * u.rad
                distribution['latitude'] = abundance.y * u.rad
                distribution['speed_dist_map'] += v_point.reshape(
                    gridlongitude.shape + (nvelbins, )) * u.s/u.km
                distribution['n_included'] += n_included.reshape(gridlongitude.shape)
                distribution['n_total'] += n_total.reshape(gridlongitude.shape)

                # Speed distribution
                velocity = mathMB.Histogram(X0_subset.speed, bins=nvelbins,
                                            range=[0, vmax], weights=weight_subset)
                distribution['speed_dist'] += velocity.histogram * u.s/u.km
                distribution['speed'] = velocity.x * u.km/u.s

                # Altitude distribution
                altitude = mathMB.Histogram(X0_subset.altitude, bins=naltbins,
                                            range=[0, np.pi / 2], weights=weight_subset)
                distribution['altitude_dist'] += altitude.histogram / u.rad
                distribution['altitude'] = altitude.x * u.rad
            
                # Azimuth distribution
                azimuth = mathMB.Histogram(X0_subset.azimuth, bins=nazbins,
                                           range=[0, 2 * np.pi], weights=weight_subset)
                distribution['azimuth_dist'] += azimuth.histogram / u.rad
                distribution['azimuth'] = azimuth.x * u.rad
                del X0, X0_subset

            distribution['fraction_observed'] = (distribution['n_included']/
                                                 distribution['n_total'])
            q = np.isnan(distribution['fraction_observed'])
            distribution['fraction_observed'][q] = 1
            distribution['abundance'] = (distribution['abundance_uncor']/
                                         distribution['fraction_observed'])
            distribution['fraction_observed'][q] = 0
            q = np.isnan(distribution['abundance'])
            distribution['abundance'][q] = 0
            
            ## normalization
            if normalize:
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

                distribution['abundance'] = (distribution['abundance'] /
                                             distribution['abundance'].sum() /
                                             area.T.to(u.cm**2) *
                                             self.sourcerate.to(1 / u.s))
                distribution['abundance_uncor'] = (distribution['abundance_uncor'] /
                                             distribution['abundance_uncor'].sum() /
                                             area.T.to(u.cm**2) *
                                             self.sourcerate.to(1 / u.s))
                
                # Normalize speed distribution
                dv = distribution['speed'][1] - distribution['speed'][0]
                distribution['speed_dist'] = (self.sourcerate *
                                              distribution['speed_dist'] /
                                              distribution['speed_dist'].sum() / dv)
                distribution['speed_dist'] = distribution['speed_dist'] * (
                    self.sourcerate.unit * u.def_unit('(km/s)^-1', u.s/u.km))
                distribution['speed_dist_map'] = distribution['speed_dist_map'] * (
                    self.sourcerate.unit * u.def_unit('(km/s)^-1', u.s/u.km))
                
                # Normalize altitude distribution
                dalt = distribution['altitude'][1] - distribution['altitude'][0]
                distribution['altitude'] = (self.sourcerate * distribution['altitude'] /
                                            distribution['altitude'].sum() / dalt)

                # Normalize azimuth distribution
                daz = distribution['azimuth'][1] - distribution['azimuth'][0]
                distribution['azimuth'] = (self.sourcerate * distribution['azimuth'] /
                                            distribution['azimuth'].sum() / daz)
            else:
                pass
                
            if which == 0:
                source = SourceMap(distribution)
                if (nlonbins > 0) and (nlatbins > 0):
                    source.abundance_uncor = distribution['abundance_uncor']
                    source.n_included = distribution['n_included']
                    source.n_total = distribution['n_total']
                    source.speed_dist_map = distribution['speed_dist_map']
                else:
                    pass
            else:
                available = SourceMap(distribution)
                if (nlonbins > 0) and (nlatbins > 0):
                    available.abundance_uncor = distribution['abundance_uncor']
                    available.n_included = distribution['n_included']
                    available.n_total = distribution['n_total']
                    available.speed_dist_map = distribution['speed_dist_map']
                else:
                    pass

        return source, available
