import os.path
import numpy as np
import pandas as pd
import pickle
import astropy.units as u
from astropy.modeling import models, fitting
from astropy.visualization import PercentileInterval
from sklearn.neighbors import KDTree
import sqlalchemy as sqla
import sqlalchemy.dialects.postgresql as pg
from nexoclom.modelcode.ModelResult import ModelResult
from nexoclom.modelcode.Output import Output
from nexoclom.modelcode.input_classes import SpatialDist, SpeedDist


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
        self.used_packets = iteration.get('used', None)
        self.used_packets0 = iteration.get('used0', None)
        

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
    
    def save(self, iteration_result, ufit_id=None):
        '''
        Insert the result of a LOS iteration into the database
        :param iteration_result: LOS result from a single outputfile
        :return: name of saved file
        '''

        # Insert the result into the database
        engine = self.inputs.config.create_engine()
        metadata_obj = sqla.MetaData()
        table = sqla.Table("uvvsmodels", metadata_obj, autoload_with=engine)
        
        insert_stmt = pg.insert(table).values(
            out_idnum = iteration_result.out_idnum,
            unfit_idnum = ufit_id,
            quantity = self.quantity,
            query = self.query,
            dphi = self.dphi,
            mechanism = self.mechanism,
            wavelength = [w.value for w in self.wavelength],
            fitted = self.fitted)
        
        with engine.connect() as con:
            result = con.execute(insert_stmt)
            con.commit()
            
        self.idnum = result.inserted_primary_key[0]
        savefile = os.path.join(os.path.dirname(iteration_result.outputfile),
                                f'model.{self.idnum}.pkl')
        print(f'Saving model result {savefile}')
        update = sqla.update(table).where(table.columns.idnum == self.idnum).values(
            filename=savefile)
        with engine.connect() as con:
            con.execute(update)
            con.commit()

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
            engine = self.inputs.config.create_engine()
            metadata_obj = sqla.MetaData()
            table = sqla.Table("uvvsmodels", metadata_obj, autoload_with=engine)
            
            query = sqla.select(table).where(
                table.columns.out_idnum == oid,
                table.columns.quantity == self.quantity,
                table.columns.query == self.query,
                table.columns.dphi == self.dphi,
                table.columns.mechanism == self.mechanism,
                table.columns.wavelength == [w.value for w in self.wavelength],
                table.columns.fitted == self.fitted)
                
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
    
    def restore(self, search_result, save_ufit_id=False):
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
    
        # Find the output files that have already been run for these inputs
        self.outid, self.outputfiles, self.npackets, self.totalsource = self.inputs.search()
        if self.npackets == 0:
            raise RuntimeError('No packets found for these Inputs.')
    
        # Find any model results that have been run for these inputs
        data = scdata.data
        search_results = self.search()
        iteration_results = []
    
        # Do this step if will need to compute any iteration results
        if None in search_results.values():
            # distance of s/c from planet
            # This is used to determine if the line of sight needs to be cut
            # short because it intersects the planet.
            dist_from_plan = np.sqrt(data.x**2 + data.y**2 + data.z**2)
        
            # Angle between look direction and planet.
            ang = np.arccos((-data.x * data.xbore - data.y * data.ybore -
                             data.z * data.zbore) / dist_from_plan)
        
            # Check to see if look direction intersects the planet anywhere
            asize_plan = np.arcsin(1. / dist_from_plan)
        
            # Don't worry about lines of sight that don't hit the planet
            dist_from_plan.loc[ang > asize_plan] = 1e30
        else:
            dist_from_plan = None
    
        for outputfile, search_result in search_results.items():
            if search_result is None:
                # simulate the data
                output = Output.restore(outputfile)
            
                packets = output.X.copy()
                packets['radvel_sun'] = (packets['vy'] +
                                         output.vrplanet.to(self.unit / u.s).value)
            
                # Note: A packet is in shadow if the line-of-sight it is on is
                #       in shadow. This is because the cone used is larger than
                #       the slit.
                self.packet_weighting(packets, output.aplanet)
            
                # This sets limits on regions where packets might be
                tree = KDTree(packets[xcols].values)
            
                rad = pd.Series(np.zeros(data.shape[0]), index=data.index)
                npack = pd.Series(np.zeros(data.shape[0]), index=data.index,
                                  dtype=int)
                used = pd.Series([set() for _ in range(data.shape[0])], index=data.index)
                used0 = pd.Series([set() for _ in range(data.shape[0])], index=data.index)
            
                print(f'{data.shape[0]} spectra taken.')
                for i, spectrum in data.iterrows():
                    x_sc = spectrum[xcols].values.astype(float)
                    bore = spectrum[borecols].values.astype(float)
                
                    # Distance from spacecraft to edge of field of view
                    a = 1
                    b = 2*np.sum(x_sc*bore)
                    c = np.linalg.norm(x_sc)**2 - self.inputs.options.outeredge**2
                    dd = (-b + np.sqrt(b**2 - 4*a*c))/2
                
                    # Compute coordinates of the LOS spaced farther apart the farther out
                    t = [np.sin(self.dphi)]
                    while t[-1] < dd:
                        t.append(t[-1] + t[-1] * np.sin(self.dphi))
                    t = np.array(t)
                    Xbore = x_sc[np.newaxis, :] + bore[np.newaxis, :] * t[:, np.newaxis]
                
                    # Narrow down number of packets
                    wid = t * np.sin(self.dphi*2)
                    ind = np.concatenate(tree.query_radius(Xbore, wid))
                    ilocs = np.unique(ind).astype(int)
                
                    subset = packets.iloc[ilocs]
                    subset_rel_sc = subset[xcols].values - x_sc[np.newaxis, :]
                    subset_dist_sc = np.linalg.norm(subset_rel_sc, axis=1)
                    losrad = np.sum(subset_rel_sc * bore[np.newaxis, :], axis=1)
                    cosang = np.sum(subset_rel_sc * bore[np.newaxis, :], axis=1)/subset_dist_sc
                    cosang[cosang > 1] = 1
                    ang = np.arccos(cosang)
                    assert np.all(np.isfinite(ang))
                
                    # Projection of packet onto line of sight
                    inview = (losrad < dist_from_plan.loc[i]) & (ang <= self.dphi)
                
                    if np.any(inview):
                        subset = subset.loc[inview]
                        subset_dist_sc = subset_dist_sc[inview]
                        losrad = losrad[inview]
                    
                        Apix = np.pi * (subset_dist_sc * np.sin(self.dphi))**2 * (
                            self.unit.to(u.cm))**2
                        wtemp = subset['weight'] / Apix
                    
                        if self.quantity == 'radiance':
                            # Determine if any packets are in shadow
                            # Projection of packet onto LOS
                            # Point along LOS the packet represents
                            hit = (x_sc[np.newaxis, :] +
                                   bore[np.newaxis, :] * losrad[:, np.newaxis])
                            rhohit = np.linalg.norm(hit[:, [0, 2]], axis=1)
                            out_of_shadow = (rhohit > 1) | (hit[:, 1] < 0)
                            wtemp *= out_of_shadow
                        
                            rad_ = wtemp.sum()
                            npack_ = np.sum(inview)
                            used.loc[i] = set(subset.loc[wtemp > 0].index)
                            used0.loc[i] = set(subset.loc[wtemp > 0, 'Index'])
                        else:
                            assert False, 'Other quantities not set up.'
                    
                        rad.loc[i] = rad_
                        npack.loc[i] = npack_
                
                    if len(data) > 10:
                        ind = data.index.get_loc(i)
                        if (ind % (len(data) // 10)) == 0:
                            print(f'Completed {ind + 1} spectra')
            
                iteration_ = {'radiance': rad,
                              'npackets': npack,
                              'totalsource': output.totalsource,
                              'outputfile': outputfile,
                              'out_idnum': output.idnum,
                              'query': scdata.query,
                              'used': used,
                              'used0': used0}
                iteration_result = IterationResult(iteration_)
                modelfile = self.save(iteration_result)
                iteration_result.modelfile = modelfile
                iteration_results.append(iteration_result)
                
                del output, packets
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
