import os.path
import numpy as np
import pandas as pd
import pickle
import random
import copy
import astropy.units as u
from sklearn.neighbors import KDTree, BallTree

from mathMB import fit_model

from .ModelResults import ModelResult
from .database_connect import database_connect
from .Input import Input
from .Output import Output


xcols = ['x', 'y', 'z']
borecols = ['xbore', 'ybore', 'zbore']
NLONBINS, NLATBINS, NVELBINS = 72, 36, 100

# Helper functions
def _should_add_weight(index, saved):
    return index in saved


def _add_weight(x, ratio):
    return np.append(x, ratio)


def _add_index(x, i):
    return np.append(x, i)


class InputError(Exception):
    """Raised when a required parameter is not included."""
    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


class LOSResult(ModelResult):
    """Class to contain the LOS result from multiple outputfiles."""
    def __init__(self, scdata, quantity='radiance', dphi=1*u.deg):
        """Determine column or emission along lines of sight.
        This assumes the model has already been run.
        
        Parameters
        ==========
        scdata
            Spacecraft data object (currently designed for MESSENGERdata object
            but can be faked for other types of data)
            
        quantity
            Quantity to calculate: 'column', 'radiance', 'density'
            
        dphi
            Angular size of the view cone. Default = r deg.
        """
        format_ = {'quantity':quantity}
        super().__init__(format_, species=scdata.species)
        
        # Basic information
        self.scdata = scdata
        self.scdata.set_frame('Model')
        self.type = 'LineOfSight'
        self.unit = None
        self.dphi = dphi.to(u.rad).value
        self.oedge = None
        
        self.fitted = None
        nspec = len(self.scdata)
        self.radiance = pd.Series(np.zeros(nspec), index=scdata.data.index)
        self.npackets = pd.Series(np.zeros(nspec), index=scdata.data.index,
                                  dtype=int)
        self.sourcemap = None
        self.modelfiles = None
        
    def delete_models(self):
        """Deletes any LOSResult models associated with this data and input
        This may never actually do anything. Overwrite=True will also
        erase the outputfiles (which erases any models that depend on them).
        Unless I put separate outputfile and modelfile delete switches,
        This shouldn't do anything"""
        
        search_results = self.search()
        if len(search_results) != 0:
            print('Warning: LOSResult.delete_models found something to delete')
            for _, search_result in search_results.items():
                if search_result is not None:
                    idnum, modelfile = search_result
                    with database_connect() as con:
                        cur = con.cursor()
                        cur.execute(f'''DELETE from uvvsmodels
                                       WHERE idnum = %s''', (idnum, ))
                    if os.path.exists(modelfile):
                        os.remove(modelfile)
        else:
            pass

    def save(self, iteration_result):
        # Insert the model into the database
        # Save is on an outputfile basis
        if self.quantity == 'radiance':
            mech = ', '.join(sorted([m for m in self.mechanism]))
            wave_ = sorted([w.value for w in self.wavelength])
            wave = ', '.join([str(w) for w in wave_])
        else:
            mech = None
            wave = None

        # Save query with all white space removed and lowercase
        tempname = f'temp_{str(random.randint(0, 1000000))}'

        with database_connect() as con:
            cur = con.cursor()
            cur.execute(f'''INSERT into uvvsmodels (out_idnum, quantity,
                            query, dphi, mechanism, wavelength,
                            fitted, filename)
                            values (%s, %s, %s, %s, %s, %s, %s, %s)''',
                        (iteration_result['out_idnum'], self.quantity,
                         self.scdata.query, self.dphi, mech, wave, self.fitted,
                         tempname))

            # Determine the savefile name
            idnum_ = pd.read_sql(f'''SELECT idnum
                                     FROM uvvsmodels
                                     WHERE filename='{tempname}';''', con)
            assert len(idnum_) == 1
            idnum = int(idnum_.idnum[0])

            savefile = os.path.join(os.path.dirname(iteration_result['outputfile']),
                                    f'model.{idnum}.pkl')

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
        for outputfile in self.outputfiles:
            with database_connect() as con:
                # Determine the id of the outputfile
                idnum_ = pd.read_sql(
                    f'''SELECT idnum
                        FROM outputfile
                        WHERE filename='{outputfile}' ''', con)
                oid = idnum_.idnum[0]
            
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
                    f'''SELECT idnum, filename FROM uvvsmodels
                        WHERE out_idnum={oid} and
                              quantity = '{self.quantity}' and
                              query = '{self.scdata.query}' and
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
                                                  result.iloc[0, 1])
                
        return search_results
    
    def restore(self, search_result):
        # Restore is on an outputfile basis
        idnum, modelfile = search_result
        with open(modelfile, 'rb') as f:
            iteration_result = pickle.load(f)
        iteration_result['idnum'] = idnum
        iteration_result['filename'] = modelfile

        return iteration_result
    
    def _data_setup(self):
        # distance of s/c from planet
        data = self.scdata.data
        dist_from_plan = np.sqrt(data.x**2 + data.y**2 + data.z**2)
    
        # Angle between look direction and planet.
        ang = np.arccos((-data.x*data.xbore - data.y*data.ybore -
                         data.z*data.zbore)/dist_from_plan)
        
        # Check to see if look direction intersects the planet anywhere
        asize_plan = np.arcsin(1./dist_from_plan)

        # Don't worry about lines of sight that don't hit the planet
        dist_from_plan.loc[ang > asize_plan] = 1e30
        
        return dist_from_plan
    
    def _spectrum_process(self, spectrum , packets, tree, dist, atoms_per_packet,
                          i=None, find_weighting=False):
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
            Apix = np.pi * (rpr[inview] * np.sin(self.dphi))**2 * (
                self.unit.to(u.cm))**2
            wtemp = subset.loc[inview, 'weight'] / Apix * atoms_per_packet
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
                pack = sum(inview)
                
                if (rad > 0) and (find_weighting):
                    ratio = spectrum.radiance / rad
    
                    # Save which packets are used for each spectrum
                    self.saved_packets[i] = subset.loc[inview, 'Index'].unique()
                    should = self.weighting.index.to_series().apply(
                        _should_add_weight, args=(self.saved_packets[i],))
                    self.weighting.loc[should] = self.weighting.loc[should].apply(
                        _add_weight, args=(ratio,))
                    self.included.loc[should] = self.included.loc[should].apply(
                        _add_index, args=(i,))
                else:
                    pass
            else:
                assert False, 'Other quantities not set up.'
        else:
            rad, pack = 0., 0
            
        return rad, pack

    def _tree(self, values, type='KDTree'):
        if type == 'KDTree':
            return KDTree(values)
        elif type == 'BallTree':
            return BallTree(values)

    def determine_source_from_data(self, inputs_, npackets, overwrite=False,
                                   packs_per_it=None,
                                   masking=None):
        self.fitted = True
        
        if isinstance(inputs_, str):
            unfit_inputs = Input(inputs_)
        elif isinstance(inputs_, Input):
            unfit_inputs = copy.deepcopy(inputs_)
        else:
            raise InputError('nexoclom.LOSResult', 'Problem with the inputs.')
        
        assert unfit_inputs.options.fitted is True, (
            'LOSResult.determine_source_from_data: inputs.options.fitted != True')

        # TAA needs to match the data
        unfit_inputs.geometry.taa = self.scdata.taa
        
        # Make the fitted inputs
        fitinput_file = os.path.join(os.path.dirname(__file__), 'data',
                                     'InputFiles',
                                     f'{self.scdata.species}.fitted.input')
        self.inputs = Input(fitinput_file)
        self.inputs.geometry.taa = self.scdata.taa
        self.unit = u.def_unit('R_' + self.inputs.geometry.planet.object,
                               self.inputs.geometry.planet.radius)

        # Run unfitted
        unfit_inputs.options.fitted = False
        unfit_inputs.run(npackets, packs_per_it=packs_per_it, overwrite=overwrite)
        unfit_outputfiles, _ , _ = unfit_inputs.search()
        
        latitude, longitude = np.ndarray((0,)), np.ndarray((0,))
        velocity, weight = np.ndarray((0,)), np.ndarray((0,))

        # Search for this result
        data = self.scdata.data
        
        iteration_results = []
        dist_from_plan = self._data_setup()
        for unfit_outputfile in unfit_outputfiles:
            # Restore the unfit output file
            output = Output.restore(unfit_outputfile)
            assert output.idnum is not None
            self.inputs.spatialdist.unfit_outid = output.idnum
            self.inputs.spatialdist.query = self.scdata.query
            self.inputs.speeddist.unfit_outid = output.idnum
            self.inputs.speeddist.query = self.scdata.query

            # Check to see if this file has already been used
            fit_outputfile, _, _ = self.inputs.search()
            assert len(fit_outputfile) <= 1
            if len(fit_outputfile) == 1:
                self.outputfiles = fit_outputfile
                search_results = self.search()
            else:
                search_results = {'default': None}
                
            for fit_outputfile, search_result in search_results.items():
                if search_result is None:
                    # need model rate for this output
                    mod_rate = output.totalsource / self.inputs.options.endtime.value
                    atoms_per_packet = 1e23 / mod_rate

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
                    # saved_packets = list of indicies of the packets used for each spectrum
                    # weighting = list of the weights that should be applied
                    #   - Final weighting for each packet is mean of weights
                    rad = pd.Series(np.zeros(data.shape[0]), index=data.index)
                    self.saved_packets = pd.Series((np.ndarray((0,), dtype=int)
                                                    for _ in range(data.shape[0])),
                                                   index=data.index)
                    ind0 = packets.Index.unique()
                    self.weighting = pd.Series((np.ndarray((0,))
                                                for _ in range(ind0.shape[0])),
                                               index=ind0)
                    self.included = pd.Series((np.ndarray((0,), dtype=np.int)
                                               for _ in range(ind0.shape[0])),
                                              index=ind0)

                    # Determine which points should be used for the fit
                    _, _, mask = fit_model(data.radiance, None, data.sigma,
                                           masking=masking, mask_only=True,
                                           altitude=data.alttan)

                    print(f'{data.shape[0]} spectra taken.')
                    for i, spectrum in data.iterrows():
                        rad_, _ = self._spectrum_process(spectrum, packets, tree,
                                                         dist_from_plan[i],
                                                         atoms_per_packet,
                                                         find_weighting=mask[i],
                                                         i=i)
                        rad.loc[i] = rad_

                        if len(data) > 10:
                            ind = data.index.get_loc(i)
                            if (ind % (len(data) // 10)) == 0:
                                print(f'Completed {ind + 1} spectra')

                    # Determine the proper weightings
                    assert np.all(self.weighting.apply(len) == self.included.apply(len))
                    new_weight = self.weighting.apply(
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
                        
                        # Save the starting state for making a source map
                        longitude = np.append(longitude, output.X0.longitude)
                        latitude = np.append(latitude, output.X0.latitude)
                        vel_ = np.sqrt(output.X0.vx**2 + output.X0.vy**2 +
                                       output.X0.vz**2) * self.inputs.geometry.planet.radius
                        velocity = np.append(velocity, vel_)
                        weight = np.append(weight, output.X0.frac)
    
                        # Update the LOSResult and output objects with new values
                        new = LOSResult(self.scdata)
                        new.fitted = True
                        new.inputs = self.inputs
                        
                        # Run the model with updated source
                        result_with_fitted = new.simulate_data_from_outputs(output)
    
                        weighting = pd.DataFrame({'weight':self.weighting.values,
                                                  'included':self.included.values})
                        iteration_result = {'radiance':result_with_fitted['radiance'],
                                            'model_total_source':output.totalsource,
                                            'atoms_per_packet':
                                                result_with_fitted['atoms_per_packet'],
                                            'weighting':weighting,
                                            'packets':self.saved_packets,
                                            'outputfile': output.filename,
                                            'out_idnum': output.idnum}
                    else:
                        iteration_result = {'radiance':pd.Series(np.zeros(data.shape[0]),
                                                                 index=data.index),
                                            'atoms_per_packet': atoms_per_packet,
                                            'model_total_source':0,
                                            'weighting':None,
                                            'packets':None,
                                            'outputfile': output.filename,
                                            'out_idnum': output.idnum}

                    print(iteration_result['outputfile'])
                    
                    modelfile = self.save(iteration_result)
                    iteration_result['modelfile'] = modelfile
                    iteration_results.append(iteration_result)
                else:
                    # Restore saved result
                    iteration_result = self.restore(search_result)
                    assert len(iteration_result['radiance']) == len(data)
                    iteration_result['modelfile'] = search_result[1]
                    iteration_results.append(iteration_result)
                    
                    output = Output.restore(fit_outputfile)
                    longitude = np.append(longitude, output.X0.longitude)
                    latitude = np.append(latitude, output.X0.latitude)
                    vel_ = np.sqrt(output.X0.vx**2 + output.X0.vy**2 +
                                   output.X0.vz**2) * self.inputs.geometry.planet.radius
                    velocity = np.append(velocity, vel_)
                    weight = np.append(weight, output.X0.frac)

        # Combine iteration_results into single new result
        self.modelfiles = {}
        for iteration_result in iteration_results:
            self.radiance += (iteration_result['radiance'] *
                              iteration_result['model_total_source'])
            self.totalsource += iteration_result['model_total_source']
            self.modelfiles[iteration_result['outputfile']] = iteration_result['modelfile']

        self.outputfiles = self.modelfiles.keys()
        self.radiance /= self.totalsource
        self.radiance *= u.R

        # Create a new sourcemap
        source, xx, yy = np.histogram2d(longitude, latitude, weights=weight,
                                        range=[[0, 2*np.pi], [-np.pi/2, np.pi/2]],
                                        bins=(NLONBINS, NLATBINS))
        source = source/np.cos(yy+(yy[1]-yy[0])/2.)[np.newaxis,:-1]
        source[:,[0,-1]] = 0
        v_source, v = np.histogram(velocity, bins=NVELBINS,
                                   range=[0, velocity.max()], weights=weight)
        v_source /= np.max(v_source)

        # packets available
        packets, _, _ = np.histogram2d(longitude, latitude,
                                       range=[[0, 2*np.pi], [-np.pi/2, np.pi/2]],
                                       bins=(NLONBINS, NLATBINS))
        packets = packets/np.cos(yy+(yy[1]-yy[0])/2.)[np.newaxis,:-1]
        packets[:,[0,-1]] = 0
        v_packets, _ = np.histogram(velocity, bins=NVELBINS, range=[0, velocity.max()])
        v_packets = v_packets / np.max(v_packets)
        
        sourcemap = {'longitude': xx*u.rad,
                     'latitude': yy*u.rad,
                     'abundance': source,
                     'p_available': packets,
                     'velocity': v*u.km/u.s,
                     'vdist': v_source,
                     'v_available': v_packets,
                     'coordinate_system': 'solar-fixed'}
        self.sourcemap = sourcemap

    def simulate_data_from_inputs(self, inputs_, npackets, overwrite=False,
                                  packs_per_it=None):
        """Given a set of inputs, determine what the spacecraft should see.
        
        Parameters
        ==========
        inputs
            A nexoclom Input object or the name of an inputs file
        """
        if isinstance(inputs_, str):
            self.inputs = Input(inputs_)
        elif isinstance(inputs_, Input):
            self.inputs = copy.deepcopy(inputs_)
        else:
            raise InputError('nexoclom.LOSResult', 'Problem with the inputs.')
        
        # TAA needs to match the data
        self.inputs.geometry.taa = self.scdata.taa
        self.unit = u.def_unit('R_' + self.inputs.geometry.planet.object,
                               self.inputs.geometry.planet.radius)

        # If using a planet-fixed source map, need to set subsolarlon
        if ((self.inputs.spatialdist.type == 'surface map') and
            (self.inputs.spatialdist.coordinate_system == 'planet-fixed')):
            self.inputs.spatialdist.subsolarlon = self.scdata.subslong.median() * u.rad
        else:
            pass
        
        # Run the model
        self.fitted = False
        self.inputs.run(npackets, packs_per_it=packs_per_it, overwrite=overwrite)
        self.search_for_outputs()

        data = self.scdata.data
        search_results = self.search()
        iteration_results = []

        dist_from_plan = (self._data_setup()
                          if None in search_results.values()
                          else None)
        for outputfile, search_result in search_results.items():
            if search_result is None:
                # simulate the data
                output = Output.restore(outputfile)
                
                # need model rate for this output
                mod_rate = output.totalsource / self.inputs.options.endtime.value
                atoms_per_packet = 1e23 / mod_rate

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
                                                         dist_from_plan[i],
                                                         atoms_per_packet)
                    rad.loc[i] = rad_
                    npack.loc[i] = pack_

                    if len(data) > 10:
                        ind = data.index.get_loc(i)
                        if (ind % (len(data)//10)) == 0:
                            print(f'Completed {ind+1} spectra')

                iteration_result = {'radiance': rad,
                                    'npackets': npack,
                                    'total_source': output.totalsource,
                                    'atoms_per_packet': atoms_per_packet,
                                    'weighting': None,
                                    'outputfile': outputfile,
                                    'out_idnum': output.idnum}
                modelfile = self.save(iteration_result)
                iteration_result['modelfile'] = modelfile
                iteration_results.append(iteration_result)
            else:
                iteration_result = self.restore(search_result)
                iteration_result['modelfile'] = search_result[1]
                assert len(iteration_result['radiance']) == len(data)
                iteration_results.append(iteration_result)

        # combine iteration_results
        self.modelfiles = {}
        for iteration_result in iteration_results:
            self.radiance += (iteration_result['radiance'] *
                              iteration_result['total_source'])
            self.npackets += iteration_result['npackets']
            self.totalsource += iteration_result['total_source']
            self.modelfiles[iteration_result['outputfile']] = iteration_result['modelfile']

        self.outputfiles = self.modelfiles.keys()
        self.radiance /= self.totalsource
        self.radiance *= u.R
        
    def simulate_data_from_outputs(self, output):
        # TAA needs to match the data
        self.inputs.geometry.taa = self.scdata.taa
        data = self.scdata.data
        self.unit = output.unit

        # need model rate for this output
        mod_rate = output.totalsource / self.inputs.options.endtime.value
        atoms_per_packet = 1e23 / mod_rate

        dist_from_plan = self._data_setup()
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
        npackets = pd.Series(np.zeros(data.shape[0]), index=data.index)
        print(f'{data.shape[0]} spectra taken.')
        for i, spectrum in data.iterrows():
            rad_, npack_ = self._spectrum_process(spectrum, packets, tree,
                                          dist_from_plan[i], atoms_per_packet)
            rad.loc[i] = rad_
            npackets.loc[i] = npack_
            
            if len(data) > 10:
                ind = data.index.get_loc(i)
                if (ind % (len(data) // 10)) == 0:
                    print(f'Completed {ind + 1} spectra')

        iteration_result = {'radiance': rad,
                            'npackets': npackets,
                            'total_source': output.totalsource,
                            'atoms_per_packet': atoms_per_packet,
                            'weighting': None,
                            'packets': None,
                            'outputfile': output.filename,
                            'out_idnum': output.idnum}
        
        return iteration_result
        
