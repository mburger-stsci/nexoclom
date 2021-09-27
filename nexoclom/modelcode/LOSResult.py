import os.path
import numpy as np
import pandas as pd
import pickle
import random
import copy
import astropy.units as u
from sklearn.neighbors import KDTree, BallTree
from nexoclom.modelcode.ModelResult import ModelResult
from nexoclom.modelcode.Output import Output
import nexoclom.math as mathMB
from nexoclom.modelcode.input_classes import SpatialDist, SpeedDist
from nexoclom.utilities import database_connect

xcols = ['x', 'y', 'z']
borecols = ['xbore', 'ybore', 'zbore']


class IterationResult:
    def __init__(self, iteration):
        self.radiance = iteration['radiance']
        self.npackets = iteration['npackets']
        self.totalsource = iteration['totalsource']
        self.outputfile = iteration['outputfile']
        self.out_idnum = iteration['out_idnum']
        self.modelfile = None
        self.model_idnum = None
        self.fitted = False
        
class FittedIterationResult(IterationResult):
    def __init__(self, iteration):
        super().__init__(iteration)
        
        self.unfit_outputfile = iteration['unfit_outputfile']
        self.unfit_outid = iteration['unfit_outid']
        self.fitted = True
        self.saved_packets = None
        self.weighting = None
        self.included = None

class InputError(Exception):
    """Raised when a required parameter is not included."""
    
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


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

        if self.fitted:
            self.unfit_outid = None
            self.unfit_outputfiles = None
        else:
            pass
        
    def __repr__(self):
        return self.__str__()
        
    def __str__(self):
        return f'''quantity = {self.quantity}
npackets = {self.npackets}
totalsource = {self.totalsource}
atoms per packet = {self.atoms_per_packet}
sourcerate = {self.sourcerate}
dphi = {self.dphi}
fit_method = {self.fit_method}
fitted = {self.fitted}'''
    
    # Helper functions
    @staticmethod
    def _should_add_weight(index, saved):
        return index in saved

    @staticmethod
    def _add_weight(x, ratio):
        return np.append(x, ratio)

    @staticmethod
    def _add_index(x, i):
        return np.append(x, i)

    # def delete_models(self):
    #     """Deletes any LOSResult models associated with this data and input
    #     This may never actually do anything. Overwrite=True will also
    #     erase the outputfiles (which erases any models that depend on them).
    #     Unless I put separate outputfile and modelfile delete switches,
    #     This shouldn't do anything"""
    #
    #     search_results = self.search()
    #     if len(search_results) != 0:
    #         print('Warning: LOSResult.delete_models found something to delete')
    #         for _, search_result in search_results.items():
    #             if search_result is not None:
    #                 idnum, modelfile = search_result
    #                 with database_connect() as con:
    #                     cur = con.cursor()
    #                     cur.execute(f'''DELETE from uvvsmodels
    #                                    WHERE idnum = %s''', (idnum,))
    #                 if os.path.exists(modelfile):
    #                     os.remove(modelfile)
    #     else:
    #         pass
    
    def save(self, iteration_result):
        '''
        Insert the result of a LOS iteration into the database
        :param iteration_result: LOS result from a single outputfile
        :return: name of saved file
        '''
        if self.quantity == 'radiance':
            mech = ', '.join(sorted([m for m in self.mechanism]))
            wave_ = sorted([w.value for w in self.wavelength])
            wave = ', '.join([str(w) for w in wave_])
        else:
            mech = None
            wave = None
        
        tempname = f'temp_{str(random.randint(0, 1000000))}'
        
        if isinstance(iteration_result, FittedIterationResult):
            ufit_id = iteration_result.unfit_outid
        else:
            ufit_id = None
        
        with database_connect() as con:
            cur = con.cursor()
            cur.execute(f'''INSERT into uvvsmodels (out_idnum, unfit_idnum,
                            quantity, query, dphi, mechanism, wavelength,
                            fitted, filename)
                            values (%s, %s, %s, %s, %s, %s, %s, %s, %s)''',
                        (iteration_result.out_idnum, ufit_id, self.quantity,
                         self.query, self.dphi, mech, wave, self.fitted,
                         tempname))
            
            # Determine the savefile name
            idnum_ = pd.read_sql(f'''SELECT idnum
                                     FROM uvvsmodels
                                     WHERE filename='{tempname}';''', con)
            assert len(idnum_) == 1
            idnum = int(idnum_.idnum[0])
            
            savefile = os.path.join(os.path.dirname(iteration_result.outputfile),
                                    f'model.{idnum}.pkl')
            print(f'Saving model result {savefile}')
            
            cur.execute(f'''UPDATE uvvsmodels
                            SET filename=%s
                            WHERE idnum=%s''', (savefile, idnum))
        
        with open(savefile, 'wb') as f:
            pickle.dump(iteration_result, f)
        
        return savefile
    
    def search(self):
        """
        :return: dictionary containing search results:
                 {outputfilename: (modelfile_id, modelfile_name)}
        """
        search_results = {}
        for oid, outputfile in zip(self.outid, self.outputfiles):
            with database_connect() as con:
                if self.quantity == 'radiance':
                    mech = ("mechanism = '" +
                            ", ".join(sorted([m for m in self.mechanism])) +
                            "'")
                    wave_ = sorted([w.value for w in self.wavelength])
                    wave = ("wavelength = '" +
                            ", ".join([str(w) for w in wave_]) +
                            "'")
                else:
                    mech = 'mechanism is NULL'
                    wave = 'wavelength is NULL'
                
                result = pd.read_sql(
                    f'''SELECT idnum, unfit_idnum, filename FROM uvvsmodels
                        WHERE out_idnum={oid} and
                              quantity = '{self.quantity}' and
                              query = '{self.query}' and
                              dphi = {self.dphi} and
                              {mech} and
                              {wave} and
                              fitted = {self.fitted}''', con)
                
                # Should only have one match per outputfile
                assert len(result) <= 1
                
                if len(result) == 0:
                    search_results[outputfile] = None
                else:
                    search_results[outputfile] = (result.iloc[0, 0],
                                                  result.iloc[0, 1],
                                                  result.iloc[0, 2])
        
        return search_results
    
    def restore(self, search_result):
        # Restore is on an outputfile basis
        idnum, ufit_idnum, modelfile = search_result
        print(f'Restoring modelfile {modelfile}.')
        with open(modelfile, 'rb') as f:
            iteration_result = pickle.load(f)
        
        iteration_result.modelfile = modelfile
        iteration_result.model_idnum = idnum
        if isinstance(iteration_result, FittedIterationResult):
            self.ufit_idnum = ufit_idnum
        else:
            pass
        
        return iteration_result
    
    def _data_setup(self, data):
        # distance of s/c from planet
        dist_from_plan = np.sqrt(data.x**2 + data.y**2 + data.z**2)
        
        # Angle between look direction and planet.
        ang = np.arccos((-data.x * data.xbore - data.y * data.ybore -
                         data.z * data.zbore) / dist_from_plan)
        
        # Check to see if look direction intersects the planet anywhere
        asize_plan = np.arcsin(1. / dist_from_plan)
        
        # Don't worry about lines of sight that don't hit the planet
        dist_from_plan.loc[ang > asize_plan] = 1e30
        
        return dist_from_plan
    
    def _spectrum_process(self, spectrum, packets, tree, dist,
                          i=None, find_weighting=False, weight_info=None):
        x_sc = spectrum[xcols].values.astype(float)
        bore = spectrum[borecols].values.astype(float)
        
        dd = 30  # Furthest distance we need to look
        x_far = x_sc + bore * dd
        while np.linalg.norm(x_far) > self.oedge:
            dd -= 0.1
            x_far = x_sc + bore * dd
        
        t = [0.05]
        while t[-1] < dd:
            t.append(t[-1] + t[-1] * np.sin(self.dphi))
        t = np.array(t)
        Xbore = x_sc[np.newaxis, :] + bore[np.newaxis, :] * t[:, np.newaxis]
        
        wid = t * np.sin(self.dphi)
        ind = np.concatenate(tree.query_radius(Xbore, wid))
        ilocs = np.unique(ind).astype(int)
        indicies = packets.iloc[ilocs].index
        subset = packets.loc[indicies]
        
        xpr = subset[xcols] - x_sc[np.newaxis, :]
        rpr = np.sqrt(xpr['x'] * xpr['x'] +
                      xpr['y'] * xpr['y'] +
                      xpr['z'] * xpr['z'])
        
        losrad = np.sum(xpr * bore[np.newaxis, :], axis=1)
        inview = rpr < dist
        
        if np.any(inview):
            # used_packets = inview.index
            # used_packets0 = packets.Index[inview.index].unique()
            
            Apix = np.pi * (rpr[inview] * np.sin(self.dphi))**2 * (
                self.unit.to(u.cm))**2
            wtemp = subset.loc[inview, 'weight'] / Apix
            if self.quantity == 'radiance':
                # Determine if any packets are in shadow
                # Projection of packet onto LOS
                # Point along LOS the packet represents
                losrad_ = losrad[inview].values
                hit = (x_sc[np.newaxis, :] +
                       bore[np.newaxis, :] * losrad_[:, np.newaxis])
                rhohit = np.linalg.norm(hit[:, [0, 2]], axis=1)
                out_of_shadow = (rhohit > 1) | (hit[:, 1] < 0)
                wtemp *= out_of_shadow
                
                rad = wtemp.sum()
                pack = np.sum(inview)
                
                if (rad > 0) and find_weighting:
                    ratio = spectrum.radiance / rad
                    
                    # Save which packets are used for each spectrum
                    weight_info['saved_packets'].loc[i] = subset.loc[inview, 'Index'].unique()
                    
                    should = weight_info['weighting'].index.to_series().apply(
                        self._should_add_weight, args=(weight_info['saved_packets'].loc[i],))
                    
                    weight_info['weighting'].loc[should] = (
                        weight_info['weighting'].loc[should].apply(self._add_weight,
                                                                   args=(ratio,)))
                    weight_info['included'].loc[should] = (
                        weight_info['included'].loc[should].apply(self._add_index,
                                                                  args=(i,)))
                else:
                    pass
            else:
                assert False, 'Other quantities not set up.'
        else:
            rad, pack, = 0., 0
            
        return rad, pack
    
    def _tree(self, values, type='KDTree'):
        if type == 'KDTree':
            return KDTree(values)
        elif type == 'BallTree':
            return BallTree(values)

    def determine_source_from_data(self, scdata):
        # Search for unfitted outputfiles
        self.fitted = True
        self.inputs.options.fitted = False
        unfit_outid, unfit_outputfiles, unfit_npackets, _ = self.inputs.search()
        if unfit_npackets == 0:
            raise RuntimeError('No packets found for these Inputs.')
        else:
            self.unfit_outid, self.unfit_outputfiles = unfit_outid, unfit_outputfiles
            self.inputs.options.fitted = True
            
            self.inputs.spatialdist = SpatialDist({'type': 'fitted output'})
            self.inputs.spatialdist.query = scdata.query
            
            self.inputs.speeddist= SpeedDist({'type': 'fitted output'})
            self.inputs.speeddist.query = scdata.query

        data = scdata.data
        iteration_results = []
        dist_from_plan = self._data_setup(data)

        # Determine which points should be used for the fit
        _, _, mask = mathMB.fit_model(data.radiance, None, data.sigma,
                                      masking=self.masking, mask_only=True,
                                      altitude=data.alttan)

        for ufit_id, ufit_out in zip(self.unfit_outid, self.unfit_outputfiles):
            # Search for fitted outputfiles
            self.inputs.spatialdist.unfit_outid = ufit_id
            self.inputs.speeddist.unfit_outid = ufit_id
            self.outid, self.outputfiles, _, _ = self.inputs.search()
            assert len(self.outid) <= 1
            
            # Search for completed fitted models
            search_results = self.search()
            if len(self.outid) == 1:
                search_result = search_results.get(self.outputfiles[0], None)
            else:
                search_result = None
            
            if search_result is None:
                output = Output.restore(ufit_out)
                
                packets = copy.deepcopy(output.X)
                packets['radvel_sun'] = (packets['vy'] +
                                         output.vrplanet.to(self.unit / u.s).value)

                self.oedge = output.inputs.options.outeredge * 2
                
                # Will base shadow on line of sight, not the packets
                out_of_shadow = np.ones(packets.shape[0])
                self.packet_weighting(packets, out_of_shadow, output.aplanet)
                
                # This sets limits on regions where packets might be
                tree = self._tree(packets[xcols].values)
                
                # rad = modeled radiance
                # weighting = list of the weights that should be applied
                #   - Final weighting for each packet is mean of weights
                rad = pd.Series(np.zeros(data.shape[0]), index=data.index)
                ind0 = packets.Index.unique()
                weight_info = {
                    'saved_packets': pd.Series((np.ndarray((0,), dtype=int)
                                                for _ in range(data.shape[0])),
                                               index=data.index),
                   'weighting': pd.Series((np.ndarray((0,))
                                           for _ in range(ind0.shape[0])),
                                          index=ind0),
                   'included': pd.Series((np.ndarray((0,), dtype=np.int)
                                          for _ in range(ind0.shape[0])),
                                          index=ind0)}
                    
                print(f'{data.shape[0]} spectra taken.')
                for i, spectrum in data.iterrows():
                    rad_, _ = self._spectrum_process(spectrum, packets, tree,
                                                     dist_from_plan[i],
                                                     find_weighting=mask[i],
                                                     i=i, weight_info=weight_info)
                    rad.loc[i] = rad_
                    
                    if len(data) > 10:
                        ind = data.index.get_loc(i)
                        if (ind % (len(data) // 10)) == 0:
                            print(f'Completed {ind + 1} spectra')
                    
                assert np.all(weight_info['weighting'].apply(len) ==
                              weight_info['included'].apply(len))

                # Determine the proper weightings
                new_weight = weight_info['weighting'].apply(
                    lambda x:x.mean() if x.shape[0] > 0 else 0.)
                new_weight /= new_weight[new_weight > 0].mean()
                assert np.all(np.isfinite(new_weight))
                
                if np.any(new_weight > 0):
                    multiplier = new_weight.loc[output.X['Index']].values
                    output.X.loc[:, 'frac'] = output.X.loc[:, 'frac'] * multiplier
                    output.X0.loc[:, 'frac'] = output.X0.loc[:, 'frac'] * new_weight
                    
                    output.X = output.X[output.X.frac > 0]
                    output.X0 = output.X0[output.X0.frac > 0]
                    output.totalsource = output.X0['frac'].sum() * output.nsteps
                    
                    # Save the fitted output
                    output.inputs = self.inputs
                    output.save()
                    
                    # Find the radiance again with the new output
                    packets = copy.deepcopy(output.X)
                    packets['radvel_sun'] = (packets['vy'] +
                                             output.vrplanet.to(self.unit / u.s).value)
                    out_of_shadow = np.ones(packets.shape[0])
                    self.packet_weighting(packets, out_of_shadow, output.aplanet)
                    tree = self._tree(packets[xcols].values)
                    rad = pd.Series(np.zeros(data.shape[0]), index=data.index)

                    for i, spectrum in data.iterrows():
                        rad_, _ = self._spectrum_process(spectrum, packets, tree,
                                                         dist_from_plan[i],
                                                         find_weighting=False,
                                                         i=i)
                        rad.loc[i] = rad_
    
                        if len(data) > 10:
                            ind = data.index.get_loc(i)
                            if (ind % (len(data) // 10)) == 0:
                                print(f'Completed {ind + 1} spectra')

                    # Save the starting state for making a source map
                    # longitude = np.append(longitude, output.X0.longitude)
                    # latitude = np.append(latitude, output.X0.latitude)
                    # vel_ = np.sqrt(output.X0.vx**2 + output.X0.vy**2 +
                    #                output.X0.vz**2) * self.inputs.geometry.planet.radius
                    # velocity = np.append(velocity, vel_)
                    # weight = np.append(weight, output.X0.frac)
                    
                    iteration = {'radiance': rad,
                                 'npackets': output.X0.frac.sum(),
                                 'totalsource': output.totalsource,
                                 'outputfile': output.filename,
                                 'out_idnum': output.idnum,
                                 'unfit_outputfile': ufit_out,
                                 'unfit_outid': ufit_id}
                    iteration_result = FittedIterationResult(iteration)
                    iteration_result.saved_packets = weight_info['saved_packets']
                    iteration_result.weighting = weight_info['weighting']
                    iteration_result.included = weight_info['included']
                else:
                    iteration = {'radiance': rad,
                                 'npackets': 0.,
                                 'totalsource': 0.,
                                 'outputfile': output.filename,
                                 'out_idnum': output.idnum,
                                 'unfit_outputfile': ufit_out,
                                 'unfit_outid': ufit_id}
                    iteration_result = FittedIterationResult(iteration)
                    iteration_result.saved_packets = weight_info['saved_packets']
                    iteration_result.weighting = weight_info['weighting']
                    iteration_result.included = weight_info['included']

                modelfile = self.save(iteration_result)
                iteration_result.modelfile = modelfile
                iteration_results.append(iteration_result)
                del output
            else:
                # Restore saved result
                print(f'Using saved file {search_result[2]}')
                iteration_result = self.restore(search_result)
                assert len(iteration_result.radiance) == len(data)
                iteration_result.model_idnum = search_result[0]
                iteration_result.modelfile = search_result[2]
                iteration_results.append(iteration_result)

        # Combine iteration_results into single new result
        self.modelfiles = {}
        for iteration_result in iteration_results:
            self.radiance += iteration_result.radiance
            self.totalsource += iteration_result.totalsource
            self.modelfiles[iteration_result.outputfile] = iteration_result.modelfile
        
        self.outputfiles = self.modelfiles.keys()
        model_rate = self.totalsource/self.inputs.options.endtime.value
        self.atoms_per_packet = 1e23 / model_rate
        self.radiance *= self.atoms_per_packet/1e3*u.kR
        self.determine_source_rate(scdata)
        self.atoms_per_packet *= self.sourcerate.unit
        self.outputfiles = list(self.modelfiles.keys())

        print(self.totalsource, self.atoms_per_packet)
    
    def simulate_data_from_inputs(self, scdata):
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
        
        # This is will work with fitted or non-fitted outputfiles
        self.outid, self.outputfiles, self.npackets, self.totalsource = self.inputs.search()
        if self.npackets == 0:
            raise RuntimeError('No packets found for these Inputs.')

        data = scdata.data
        search_results = self.search()
        iteration_results = []
        
        # Do this step if will need to compute any iteration results
        dist_from_plan = (self._data_setup(data)
                          if None in search_results.values()
                          else None)
        for outputfile, search_result in search_results.items():
            if search_result is None:
                # simulate the data
                output = Output.restore(outputfile)
                
                packets = copy.deepcopy(output.X)
                packets['radvel_sun'] = (packets['vy'] +
                                         output.vrplanet.to(self.unit / u.s).value)
                self.oedge = output.inputs.options.outeredge * 2
                
                # Will base shadow on line of sight, not the packets
                out_of_shadow = np.ones(len(packets))
                self.packet_weighting(packets, out_of_shadow, output.aplanet)
                
                # This sets limits on regions where packets might be
                tree = self._tree(packets[xcols].values)
                
                rad = pd.Series(np.zeros(data.shape[0]), index=data.index)
                npack = pd.Series(np.zeros(data.shape[0]), index=data.index,
                                  dtype=int)
                print(f'{data.shape[0]} spectra taken.')
                for i, spectrum in data.iterrows():
                    rad_, pack_ = self._spectrum_process(spectrum, packets, tree,
                                                         dist_from_plan[i])
                    rad.loc[i] = rad_
                    npack.loc[i] = pack_
                    
                    if len(data) > 10:
                        ind = data.index.get_loc(i)
                        if (ind % (len(data) // 10)) == 0:
                            print(f'Completed {ind + 1} spectra')
                
                iteration_ = {'radiance': rad,
                              'npackets': npack,
                              'totalsource': output.totalsource,
                              'outputfile': outputfile,
                              'out_idnum': output.idnum,
                              'query': scdata.query}
                iteration_result = IterationResult(iteration_)
                modelfile = self.save(iteration_result)
                iteration_result.modelfile = modelfile
                iteration_results.append(iteration_result)
            else:
                print(f'Using saved result {search_result[2]}')
                iteration_result = self.restore(search_result)
                iteration_result.model_idnum = search_result[0]
                iteration_result.modelfile = search_result[2]
                assert len(iteration_result.radiance) == len(data)
                iteration_results.append(iteration_result)

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
        self.atoms_per_packet *= self.sourcerate.unit
        self.outputfiles = self.modelfiles.keys()
        
        print(self.totalsource, self.atoms_per_packet)
        
        # self.sourcemap = self.make_source_map(longitude, latitude, velocity,
        #                                       weight)
        
    def determine_source_rate(self, scdata):
        strength, goodness_of_fit, mask = mathMB.fit_model(scdata.data.radiance.values,
                                                           self.radiance.values,
                                                           scdata.data.sigma.values,
                                                           fit_method=self.fit_method,
                                                           masking=self.masking,
                                                           altitude=scdata.data.alttan)
        self.radiance *= strength
        self.sourcerate = strength * u.def_unit('10**23 atoms/s', 1e23 / u.s)
        self.goodness_of_fit = goodness_of_fit
        self.mask = mask
        